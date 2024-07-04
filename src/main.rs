mod hittables;

use std::f32::consts::PI;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufWriter;
use std::ops::Range;

use hittables::Hittable;
use hittables::HittableList;
use raylib::math::Vector3;

const EPS: f32 = 1e-3;

const FILENAME: &str = "image.ppm";
const SCREEN_FACTOR: u32 = 200;
const IMAGE_WIDTH: u32 = 3 * SCREEN_FACTOR;
const IMAGE_HEIGHT: u32 = 2 * SCREEN_FACTOR;
const COLOR_MAX: f32 = 256.0;

const VFOV: f32 = 20.0;

const LOOK_FROM: Point3 = Point3::new(13.0, 2.0, 3.0);
const LOOK_AT: Point3 = Point3::new(0.0, 0.0, 0.0);
const CAMERA_UP: Vector3 = Vector3::new(0.0, 1.0, 0.0);

const DEFOCUS_ANGLE: f32 = 0.0;
const FOCUS_DIST: f32 = 1.0;

const SAMPLES_PER_PIXEL: u32 = 10;
const MAX_SAMPLE_DEPTH: u32 = 50;

const PROGRESS_BAR_SEGMENTS: u32 = 20;

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

fn ray_color(ray: &Ray, hittable: &impl Hittable, depth: u32) -> Color {
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

#[derive(Clone, Copy)]
pub enum Material {
    Lambertian { albedo: Color },
    Metal { albedo: Color, fuzz: f32 },
    Dielectric { refraction_index: f32 },
}

impl Material {
    const fn lambertian(albedo: Color) -> Material {
        Material::Lambertian { albedo }
    }
    const fn metal(albedo: Color, fuzz: f32) -> Material {
        Material::Metal { albedo, fuzz }
    }
    const fn dielectric(refraction_index: f32) -> Material {
        Material::Dielectric { refraction_index }
    }

    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Color, Ray)> {
        match self {
            Self::Lambertian { albedo } => self.scatter_lambertian(ray, hit, *albedo),
            Self::Metal { albedo, fuzz } => self.scatter_metal(ray, hit, *albedo, *fuzz),
            Self::Dielectric { refraction_index } => {
                self.scatter_dielectric(ray, hit, *refraction_index)
            }
        }
    }

    fn scatter_lambertian(&self, _ray: &Ray, hit: &Hit, albedo: Color) -> Option<(Color, Ray)> {
        let mut scatter_dir = hit.normal + random_unit_vector();
        if vec3_near_zero(scatter_dir) {
            scatter_dir = hit.normal
        }
        Some((albedo, Ray::new(hit.p, scatter_dir)))
    }

    fn scatter_metal(
        &self,
        ray: &Ray,
        hit: &Hit,
        albedo: Color,
        fuzz: f32,
    ) -> Option<(raylib::prelude::Vector3, Ray)> {
        let mut reflected = reflect(ray.dir, hit.normal);
        reflected = reflected.normalized() + random_unit_vector() * fuzz;

        // If it scattered backwards, don't
        if hit.normal.dot(reflected) <= 0.0 {
            return None;
        }
        Some((albedo, Ray::new(hit.p, reflected)))
    }

    fn scatter_dielectric(
        &self,
        ray: &Ray,
        hit: &Hit,
        refraction_index: f32,
    ) -> Option<(Color, Ray)> {
        let ri = if hit.front_face {
            1.0 / refraction_index
        } else {
            refraction_index
        };

        let unit_dir = ray.dir.normalized();
        let cos_theta = -unit_dir.dot(hit.normal).min(1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        let cannot_refract = ri * sin_theta > 1.0;

        let direction = if cannot_refract || reflectance(cos_theta, ri) > fastrand::f32() {
            reflect(unit_dir, hit.normal)
        } else {
            refract(unit_dir, hit.normal, ri)
        };

        Some((Color::one(), Ray::new(hit.p, direction)))
    }
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

fn vec3_near_zero(v: Vector3) -> bool {
    let thresh = 1e-8;
    v.x.abs() < thresh && v.y.abs() < thresh && v.z.abs() < thresh
}

fn random_f32(min: f32, max: f32) -> f32 {
    fastrand::f32() * (max - min) + min
}

fn random_vec3() -> Vector3 {
    Vector3::new(fastrand::f32(), fastrand::f32(), fastrand::f32())
}

fn random_vec3_range(min: f32, max: f32) -> Vector3 {
    random_vec3() * (max - min) + min
}

fn random_vec3_in_sphere(radius: f32) -> Vector3 {
    loop {
        let candidate = random_vec3_range(-radius, radius);
        if candidate.dot(candidate) < radius * radius {
            return candidate;
        }
    }
}

fn random_unit_vector() -> Vector3 {
    random_vec3_in_sphere(1.0).normalized()
}

fn random_vector_in_unit_disk() -> Vector3 {
    loop {
        let mut candidate = random_vec3_range(-1.0, 1.0);
        candidate.z = 0.0;

        if candidate.dot(candidate) < 1.0 {
            return candidate;
        }
    }
}

fn defocus_disk_sample(center: Point3, defocus_disk_u: Vector3, defocus_disk_v: Vector3) -> Point3 {
    let p = random_vector_in_unit_disk();
    center + defocus_disk_u * p.x + defocus_disk_v * p.y
}

fn reflect(v: Vector3, n: Vector3) -> Vector3 {
    v - n * v.dot(n) * 2.0
}

fn refract(uv: Vector3, n: Vector3, eta_i_over_eta_t: f32) -> Vector3 {
    let cos_theta = -uv.dot(n).min(1.0);
    let r_out_perp = (uv + n * cos_theta) * eta_i_over_eta_t;
    let r_out_parallel = n * -((1.0 - r_out_perp.dot(r_out_perp)).abs()).sqrt();
    r_out_perp + r_out_parallel
}

fn reflectance(cosine: f32, refraction_index: f32) -> f32 {
    let mut r0 = (1.0 - refraction_index) / (1.0 + refraction_index);
    r0 *= r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
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
    let world = make_world();

    // Write output image file
    let f = File::create(FILENAME)?;
    {
        let mut w = BufWriter::new(f);

        // Write file header
        writeln!(&mut w, "P3")?;
        writeln!(&mut w, "{} {}", IMAGE_WIDTH, IMAGE_HEIGHT)?;
        writeln!(&mut w, "{}", COLOR_MAX as u32 - 1)?;

        let mut progress_bar = "".to_string();
        println!("Processing image...");

        for j in 0..IMAGE_HEIGHT {
            if (j + 1) % (IMAGE_HEIGHT / PROGRESS_BAR_SEGMENTS) == 0 {
                progress_bar += "=";
            }
            print!(
                "\r\x1b[K[{: <1$}]",
                progress_bar, PROGRESS_BAR_SEGMENTS as usize
            );
            print!(" {}/{}", j + 1, IMAGE_HEIGHT);
            io::stdout().flush().unwrap();
            for i in 0..IMAGE_WIDTH {
                let pixel_center = p0 + du * i as f32 + dv * j as f32;

                let mut c = Color::zero();
                for _ in 0..SAMPLES_PER_PIXEL {
                    let offset = Vector3::new(fastrand::f32() - 0.5, fastrand::f32() - 0.5, 0.0);
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
                write_color(&mut w, c)?;
            }
        }
        println!("\nDone!")
    }

    Ok(())
}
