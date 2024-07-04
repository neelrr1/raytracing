mod hittables;
mod material;
mod utils;

use std::f32::consts::PI;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufWriter;
use std::ops::Range;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::thread;

use hittables::HittableList;
use material::Material;
use raylib::math::Vector3;
use utils::defocus_disk_sample;
use utils::random_f32;
use utils::random_vec3;
use utils::random_vec3_range;

const EPS: f32 = 1e-3;

const FILENAME: &str = "image.ppm";
const SCREEN_FACTOR: usize = 400;
const IMAGE_WIDTH: usize = 3 * SCREEN_FACTOR;
const IMAGE_HEIGHT: usize = 2 * SCREEN_FACTOR;
const COLOR_MAX: f32 = 256.0;

const VFOV: f32 = 20.0;

const LOOK_FROM: Point3 = Point3::new(13.0, 2.0, 3.0);
const LOOK_AT: Point3 = Point3::new(0.0, 0.0, 0.0);
const CAMERA_UP: Vector3 = Vector3::new(0.0, 1.0, 0.0);

const DEFOCUS_ANGLE: f32 = 0.2;
const FOCUS_DIST: f32 = 1.0;

const SAMPLES_PER_PIXEL: u32 = 500;
const MAX_SAMPLE_DEPTH: u32 = 50;

const NUM_THREADS: usize = 10;

type Color = Vector3;
type Point3 = Vector3;

// Gamma correction for colors
fn linear_to_gamma(x: f32) -> f32 {
    if x > 0.0 {
        x.sqrt()
    } else {
        0.0
    }
}

fn write_color(mut w: impl Write, mut c: Color) -> io::Result<()> {
    // Gamma correction
    c.x = linear_to_gamma(c.x);
    c.y = linear_to_gamma(c.y);
    c.z = linear_to_gamma(c.z);

    c *= COLOR_MAX - EPS;

    writeln!(&mut w, "{} {} {}", c.x as u32, c.y as u32, c.z as u32)?;
    Ok(())
}

fn write_image_file(pixels: Vec<Color>) -> io::Result<()> {
    let f = File::create(FILENAME)?;
    {
        let mut w = BufWriter::new(f);

        // Write file header
        writeln!(&mut w, "P3")?;
        writeln!(&mut w, "{} {}", IMAGE_WIDTH, IMAGE_HEIGHT)?;
        writeln!(&mut w, "{}", COLOR_MAX as u32 - 1)?;

        for pixel in pixels {
            write_color(&mut w, pixel)?;
        }
    }
    Ok(())
}

pub struct Ray {
    origin: Point3,
    dir: Point3,
}

impl Ray {
    fn new(origin: Point3, dir: Point3) -> Ray {
        Ray { origin, dir }
    }

    fn at(&self, t: f32) -> Point3 {
        self.origin + self.dir * t
    }
}

fn ray_color(ray: &Ray, hittable: &HittableList, depth: u32) -> Color {
    if depth == 0 {
        return Color::zero();
    }

    if let Some(hit) = hittable.hit(ray, EPS, f32::INFINITY) {
        if let Some((attenuation, scattered)) = hit.mat.scatter(ray, &hit) {
            return attenuation * ray_color(&scattered, hittable, depth - 1);
        }
        return Color::zero();
    }

    let d = ray.dir.normalized();

    // remap y from 0 to 1
    let a = (d.y + 1.0) * 0.5;
    Color::one().lerp(Color::new(0.5, 0.7, 1.0), a)
}

pub struct Hit {
    p: Point3,
    normal: Vector3,
    t: f32,
    mat: Material,
    front_face: bool,
}

impl Hit {
    fn new(p: Point3, normal: Vector3, t: f32, mat: Material, front_face: bool) -> Hit {
        Hit {
            p,
            normal,
            t,
            mat,
            front_face,
        }
    }
}

const GROUND: Material = Material::lambertian(Color::new(0.5, 0.5, 0.5));
const GLASS: Material = Material::dielectric(1.5);

fn make_world() -> HittableList {
    let mut world = HittableList::new();

    world.add_sphere(Point3::new(0.0, -1000.0, 0.0), 1000.0, GROUND);
    for i in -11..11 {
        for j in -11..11 {
            let mat_choice = fastrand::f32();
            let center = Point3::new(
                i as f32 + 0.9 * fastrand::f32(),
                0.2,
                j as f32 + 0.9 * fastrand::f32(),
            );

            if (center - Point3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                if mat_choice < 0.8 {
                    let mat = Material::lambertian(random_vec3() * random_vec3());
                    world.add_sphere(center, 0.2, mat);
                } else if mat_choice < 0.95 {
                    let mat = Material::metal(random_vec3_range(0.5, 1.0), random_f32(0.0, 0.5));
                    world.add_sphere(center, 0.2, mat);
                } else {
                    world.add_sphere(center, 0.2, GLASS);
                }
            }
        }
    }
    world.add_sphere(Point3::up(), 1.0, GLASS);

    let mat2 = Material::lambertian(Color::new(0.4, 0.2, 0.1));
    world.add_sphere(Point3::new(-4.0, 1.0, 0.0), 1.0, mat2);

    let mat3 = Material::metal(Color::new(0.7, 0.6, 0.5), 0.0);
    world.add_sphere(Point3::new(4.0, 1.0, 0.0), 1.0, mat3);

    world
}

fn main() -> io::Result<()> {
    // Calculate viewport dimensions
    let theta = VFOV / 180.0 * PI;
    let h = (theta / 2.0).tan();
    let viewport_height = 2.0 * h * FOCUS_DIST;
    let viewport_width = IMAGE_WIDTH as f32 / IMAGE_HEIGHT as f32 * viewport_height;

    // Calculate unit vectors relative to camera
    let w = (LOOK_FROM - LOOK_AT).normalized();
    let u = CAMERA_UP.cross(w).normalized();
    let v = w.cross(u);

    let viewport_u = u * viewport_width;
    let viewport_v = -v * viewport_height;

    // Camera
    let camera_center = LOOK_FROM;
    let viewport_upper_left: Vector3 =
        camera_center - w * FOCUS_DIST - viewport_u / 2.0 - viewport_v / 2.0;
    let du = viewport_u / IMAGE_WIDTH as f32;
    let dv = viewport_v / IMAGE_HEIGHT as f32;
    let p0 = viewport_upper_left + du / 2.0 + dv / 2.0;

    // Basis vectors for defocus disk
    let defocus_angle = (DEFOCUS_ANGLE / 2.0) * PI / 180.0;
    let defocus_radius = FOCUS_DIST * defocus_angle.tan();
    let defocus_disk_u = u * defocus_radius;
    let defocus_disk_v = v * defocus_radius;

    // Set up the scene
    let world = Arc::new(make_world());

    // Perform raytracing
    println!("Processing...");

    let mut children = Vec::new();
    let mut chunk_start = 0;
    let rem = IMAGE_HEIGHT % NUM_THREADS;
    let mut chunk_size = IMAGE_HEIGHT;

    let finished = Arc::new(AtomicUsize::new(0));

    for i in 0..NUM_THREADS {
        let world = Arc::clone(&world);
        let finished = Arc::clone(&finished);

        if NUM_THREADS > 1 {
            chunk_size = IMAGE_HEIGHT / NUM_THREADS;
            if i < rem {
                chunk_size += 1;
            }
        }

        children.push(thread::spawn(move || {
            let mut pixels = Vec::new();

            for j in chunk_start..(chunk_start + chunk_size) {
                print!(
                    "\r\x1b[K{}/{}",
                    finished.load(Ordering::Relaxed) + 1,
                    IMAGE_HEIGHT
                );
                io::stdout().flush().unwrap();
                for i in 0..IMAGE_WIDTH {
                    let pixel_center = p0 + du * i as f32 + dv * j as f32;

                    let mut c = Color::zero();
                    for _ in 0..SAMPLES_PER_PIXEL {
                        let offset =
                            Vector3::new(fastrand::f32() - 0.5, fastrand::f32() - 0.5, 0.0);
                        let sample_center = pixel_center + offset * du + offset * dv;

                        // not normalized
                        let dir = sample_center - camera_center;
                        let r = Ray::new(
                            if defocus_angle <= 0.0 {
                                camera_center
                            } else {
                                defocus_disk_sample(camera_center, defocus_disk_u, defocus_disk_v)
                            },
                            dir,
                        );

                        c += ray_color(&r, &world, MAX_SAMPLE_DEPTH);
                    }

                    c /= SAMPLES_PER_PIXEL as f32;
                    c.clamp(Range {
                        start: 0.0,
                        end: COLOR_MAX,
                    });
                    pixels.push(c);
                }
                finished.fetch_add(1, Ordering::Relaxed);
            }
            pixels
        }));
        chunk_start += chunk_size;
    }

    // Collect results from all threads
    let pixels =
        children
            .into_iter()
            .map(|c| c.join().unwrap())
            .fold(Vec::new(), |mut acc, mut e| {
                acc.append(&mut e);
                acc
            });
    assert!(pixels.len() == IMAGE_HEIGHT * IMAGE_WIDTH);

    // Write output image file
    write_image_file(pixels)?;

    println!();
    println!();
    println!("Done!");

    Ok(())
}
