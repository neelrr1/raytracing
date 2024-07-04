use crate::{
    utils::{random_unit_vector, reflect, reflectance, refract, vec3_near_zero},
    Color, Hit, Ray,
};

#[derive(Clone, Copy)]
pub enum Material {
    Lambertian { albedo: Color },
    Metal { albedo: Color, fuzz: f32 },
    Dielectric { refraction_index: f32 },
}

impl Material {
    pub const fn lambertian(albedo: Color) -> Material {
        Material::Lambertian { albedo }
    }
    pub const fn metal(albedo: Color, fuzz: f32) -> Material {
        Material::Metal { albedo, fuzz }
    }
    pub const fn dielectric(refraction_index: f32) -> Material {
        Material::Dielectric { refraction_index }
    }

    pub fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Color, Ray)> {
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
