use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub struct Vector3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vector3 {
    pub const fn new(x: f32, y: f32, z: f32) -> Vector3 {
        Vector3 { x, y, z }
    }

    /// Returns a new `Vector3` with all components set to zero.
    pub fn zero() -> Vector3 {
        Vector3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Returns a new `Vector3` with all components set to one.
    pub fn one() -> Vector3 {
        Vector3 {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        }
    }

    /// Returns a new `Vector3` containing the cross product between `self` and vector `v`.
    pub fn cross(&self, v: &Vector3) -> Vector3 {
        Vector3 {
            x: self.y * v.z - self.z * v.y,
            y: self.z * v.x - self.x * v.z,
            z: self.x * v.y - self.y * v.x,
        }
    }

    /// Returns a new `Vector3` with componenets linearly interpolated by `amount` towards vector `v`.
    pub fn lerp(&self, v: &Vector3, amount: f32) -> Vector3 {
        Vector3 {
            x: self.x + amount * (v.x - self.x),
            y: self.y + amount * (v.y - self.y),
            z: self.z + amount * (v.z - self.z),
        }
    }

    /// Returns a new `Vector3` with componenets clamp to a certain interval.
    pub fn clamp(&self, start: f32, end: f32) -> Vector3 {
        Vector3 {
            x: self.x.clamp(start, end),
            y: self.y.clamp(start, end),
            z: self.z.clamp(start, end),
        }
    }

    /// Calculates the vector length.
    pub fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Calculates the dot product with vector `v`.
    pub fn dot(&self, v: &Vector3) -> f32 {
        self.x * v.x + self.y * v.y + self.z * v.z
    }

    /// Returns a new `Vector3` with normalized components from the current vector.
    pub fn normalized(&self) -> Vector3 {
        let mut length = self.length();
        if length == 0.0 {
            length = 1.0;
        }
        let ilength = 1.0 / length;

        Vector3 {
            x: self.x * ilength,
            y: self.y * ilength,
            z: self.z * ilength,
        }
    }
}

impl Add for Vector3 {
    type Output = Vector3;
    fn add(self, v: Vector3) -> Self {
        Vector3 {
            x: self.x + v.x,
            y: self.y + v.y,
            z: self.z + v.z,
        }
    }
}

impl Add<f32> for Vector3 {
    type Output = Vector3;
    fn add(self, value: f32) -> Self {
        Vector3 {
            x: self.x + value,
            y: self.y + value,
            z: self.z + value,
        }
    }
}

impl AddAssign for Vector3 {
    fn add_assign(&mut self, v: Vector3) {
        self.x += v.x;
        self.y += v.y;
        self.z += v.z;
    }
}

impl AddAssign<f32> for Vector3 {
    fn add_assign(&mut self, value: f32) {
        self.x += value;
        self.y += value;
        self.z += value;
    }
}

impl Sub for Vector3 {
    type Output = Vector3;
    fn sub(self, v: Vector3) -> Self {
        Vector3 {
            x: self.x - v.x,
            y: self.y - v.y,
            z: self.z - v.z,
        }
    }
}

impl Sub<f32> for Vector3 {
    type Output = Vector3;
    fn sub(self, value: f32) -> Self {
        Vector3 {
            x: self.x - value,
            y: self.y - value,
            z: self.z - value,
        }
    }
}

impl SubAssign for Vector3 {
    fn sub_assign(&mut self, v: Vector3) {
        self.x -= v.x;
        self.y -= v.y;
        self.z -= v.z;
    }
}

impl SubAssign<f32> for Vector3 {
    fn sub_assign(&mut self, value: f32) {
        self.x -= value;
        self.y -= value;
        self.z -= value;
    }
}

impl Mul for Vector3 {
    type Output = Vector3;
    fn mul(self, v: Vector3) -> Self {
        Vector3 {
            x: self.x * v.x,
            y: self.y * v.y,
            z: self.z * v.z,
        }
    }
}

impl Mul<f32> for Vector3 {
    type Output = Vector3;
    fn mul(self, value: f32) -> Self {
        Vector3 {
            x: self.x * value,
            y: self.y * value,
            z: self.z * value,
        }
    }
}

impl MulAssign for Vector3 {
    fn mul_assign(&mut self, v: Vector3) {
        self.x *= v.x;
        self.y *= v.y;
        self.z *= v.z;
    }
}

impl MulAssign<f32> for Vector3 {
    fn mul_assign(&mut self, value: f32) {
        self.x *= value;
        self.y *= value;
        self.z *= value;
    }
}

impl Div for Vector3 {
    type Output = Vector3;
    fn div(self, v: Vector3) -> Self {
        Vector3 {
            x: self.x / v.x,
            y: self.y / v.y,
            z: self.z / v.z,
        }
    }
}

impl Div<f32> for Vector3 {
    type Output = Vector3;
    fn div(self, value: f32) -> Self {
        Vector3 {
            x: self.x / value,
            y: self.y / value,
            z: self.z / value,
        }
    }
}

impl DivAssign for Vector3 {
    fn div_assign(&mut self, v: Vector3) {
        self.x /= v.x;
        self.y /= v.y;
        self.z /= v.z;
    }
}

impl DivAssign<f32> for Vector3 {
    fn div_assign(&mut self, value: f32) {
        self.x /= value;
        self.y /= value;
        self.z /= value;
    }
}

impl Neg for Vector3 {
    type Output = Vector3;
    fn neg(self) -> Self {
        Vector3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
