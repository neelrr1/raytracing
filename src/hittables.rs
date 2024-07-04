use crate::{Hit, Material, Point3, Ray};

pub enum Hittable {
    Sphere {
        center: Point3,
        radius: f32,
        mat: Material,
    },
}

impl Hittable {
    pub fn sphere(center: Point3, radius: f32, mat: Material) -> Hittable {
        Hittable::Sphere {
            center,
            radius,
            mat,
        }
    }

    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
        match self {
            Hittable::Sphere {
                center,
                radius,
                mat,
            } => Hittable::hit_sphere(ray, t_min, t_max, center, radius, mat),
        }
    }

    fn hit_sphere(
        ray: &Ray,
        t_min: f32,
        t_max: f32,
        center: &Point3,
        radius: &f32,
        mat: &Material,
    ) -> Option<Hit> {
        let ray_to_sphere = *center - ray.origin;
        let a = ray.dir.dot(&ray.dir);
        let h = ray.dir.dot(&ray_to_sphere);
        let c = ray_to_sphere.dot(&ray_to_sphere) - radius * radius;

        let discrim = h * h - a * c;

        let sqrtd = discrim.sqrt();
        let mut root = (h - sqrtd) / a;
        if root.is_nan() || root < t_min || root > t_max {
            // Check if the other root is in range
            root = (h + sqrtd) / a;
            if root < t_min || root > t_max || root.is_nan() {
                return None;
            }
        }

        let rp = ray.at(root);
        let normal = (rp - *center) / *radius;
        let front_face = ray.dir.dot(&normal) < 0.0;

        let out = Hit::new(
            rp,
            if front_face { normal } else { -normal },
            root,
            *mat,
            front_face,
        );
        Some(out)
    }
}

pub struct HittableList {
    objects: Vec<Hittable>,
}

impl HittableList {
    pub fn new() -> HittableList {
        HittableList {
            objects: Vec::new(),
        }
    }

    pub fn add_object(&mut self, o: Hittable) {
        self.objects.push(o);
    }

    pub fn add_sphere(&mut self, center: Point3, radius: f32, mat: Material) {
        self.add_object(Hittable::sphere(center, radius, mat))
    }

    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<Hit> {
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
