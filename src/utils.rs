use crate::{vector3::Vector3, Point3};

pub fn vec3_near_zero(v: Vector3) -> bool {
    let thresh = 1e-8;
    v.x.abs() < thresh && v.y.abs() < thresh && v.z.abs() < thresh
}

pub fn random_f32(min: f32, max: f32) -> f32 {
    fastrand::f32() * (max - min) + min
}

pub fn random_vec3() -> Vector3 {
    Vector3::new(fastrand::f32(), fastrand::f32(), fastrand::f32())
}

pub fn random_vec3_range(min: f32, max: f32) -> Vector3 {
    random_vec3() * (max - min) + min
}

pub fn random_vec3_in_sphere(radius: f32) -> Vector3 {
    loop {
        let candidate = random_vec3_range(-radius, radius);
        if candidate.dot(&candidate) < radius * radius {
            return candidate;
        }
    }
}

pub fn random_unit_vector() -> Vector3 {
    random_vec3_in_sphere(1.0).normalized()
}

pub fn random_vector_in_unit_disk() -> Vector3 {
    loop {
        let mut candidate = random_vec3_range(-1.0, 1.0);
        candidate.z = 0.0;

        if candidate.dot(&candidate) < 1.0 {
            return candidate;
        }
    }
}

pub fn defocus_disk_sample(
    center: Point3,
    defocus_disk_u: Vector3,
    defocus_disk_v: Vector3,
) -> Point3 {
    let p = random_vector_in_unit_disk();
    center + defocus_disk_u * p.x + defocus_disk_v * p.y
}

pub fn reflect(v: Vector3, n: Vector3) -> Vector3 {
    v - n * v.dot(&n) * 2.0
}

pub fn refract(uv: Vector3, n: Vector3, eta_i_over_eta_t: f32) -> Vector3 {
    let cos_theta = -uv.dot(&n).min(1.0);
    let r_out_perp = (uv + n * cos_theta) * eta_i_over_eta_t;
    let r_out_parallel = n * -((1.0 - r_out_perp.dot(&r_out_perp)).abs()).sqrt();
    r_out_perp + r_out_parallel
}

pub fn reflectance(cosine: f32, refraction_index: f32) -> f32 {
    let mut r0 = (1.0 - refraction_index) / (1.0 + refraction_index);
    r0 *= r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}
