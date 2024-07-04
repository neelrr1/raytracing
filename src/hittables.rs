use crate::{Hit, Material, Point3, Ray};

pub trait Hittable {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<Hit>;
}

pub struct Sphere {
    center: Point3,
    radius: f32,
    mat: &'static dyn Material,
}

impl Sphere {
    pub fn new(center: Point3, radius: f32, mat: &'static impl Material) -> Sphere {
        Sphere {
            center,
            radius,
            mat,
        }
    }
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
        let ray_to_sphere = self.center - ray.origin;
        let a = ray.dir.dot(ray.dir);
        let h = ray.dir.dot(ray_to_sphere);
        let c = ray_to_sphere.dot(ray_to_sphere) - self.radius * self.radius;

        let discrim = h * h - a * c;

        let sqrtd = discrim.sqrt();
        let mut root = (h - sqrtd) / a;
        if root.is_nan() || root < t_min || root > t_max {
            // Check if the other root is in range
            root = (h + sqrtd) / a;
            if root.is_nan() || root < t_min || root > t_max {
                return None;
            }
        }

        let rp = ray.at(root);
        let normal = (rp - self.center) / self.radius;
        let front_face = ray.dir.dot(normal) < 0.0;

        let out = Hit::new(
            rp,
            if front_face { normal } else { -normal },
            root,
            self.mat,
            front_face,
        );
        Some(out)
    }
}

pub struct HittableList {
    objects: Vec<Box<dyn Hittable>>,
}

impl HittableList {
    pub fn new() -> HittableList {
        HittableList {
            objects: Vec::new(),
        }
    }

    pub fn add_object(&mut self, o: impl Hittable + 'static) {
        self.objects.push(Box::new(o));
    }

    pub fn add_sphere(&mut self, center: Point3, radius: f32, mat: &'static impl Material) {
        self.add_object(Sphere::new(center, radius, mat))
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
        let mut hit = None;
        let mut t = t_max;

        for hittable in &self.objects {
            if let Some(new_hit) = hittable.hit(ray, t_min, t_max) {
                if new_hit.t < t {
                    t = new_hit.t;
                    hit = Some(new_hit);
                }
            }
        }

        hit
    }
}
