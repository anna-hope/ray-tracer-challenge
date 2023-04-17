use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::equal;

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        equal(self.x, other.x) && equal(self.y, other.y) && equal(self.z, other.z)
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Add<Vector> for Point {
    type Output = Self;

    fn add(self, rhs: Vector) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Mul<f64> for Point {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<f64> for Point {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Sub for Point {
    type Output = Vector;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Sub<Vector> for Point {
    type Output = Self;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn magnitude(&self) -> f64 {
        let sum = self.x.powi(2) + self.y.powi(2) + self.z.powi(2);
        sum.sqrt()
    }

    pub fn norm(&self) -> Self {
        let magnitude = self.magnitude();
        Self {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
        }
    }

    pub fn dot(&self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn reflect(self, normal: Self) -> Self {
        self - normal * 2. * self.dot(normal)
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        equal(self.x, other.x) && equal(self.y, other.y) && equal(self.z, other.z)
    }
}

impl Add for Vector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Add<Point> for Vector {
    type Output = Self;

    fn add(self, rhs: Point) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f64> for Vector {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<f64> for Vector {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Neg for Vector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn point_makes_point() {
        let _a = Point::new(4.0, -4.0, 3.0);
    }

    #[test]
    fn vector_makes_vector() {
        let _a = Vector::new(4.0, -4.0, 3.0);
    }

    #[test]
    fn points_and_vectors_equal() {
        let a = Point::new(1.0, 1.0, 1.0);
        let b = Point::new(1.0, 1.0, 1.0);
        assert_eq!(a, b);

        let a = Vector::new(1.0, 2.0, -3.0);
        let b = Vector::new(1.0, 2.0, -3.0);
        assert_eq!(a, b);
    }

    #[test]
    fn points_and_vectors_not_equal() {
        let a = Point::new(1.0, 1.0, -1.0);
        let b = Point::new(1.0, 1.0, 1.0);
        assert_ne!(a, b);

        let a = Vector::new(1.0, 2.0, 3.0);
        let b = Vector::new(1.0, 2.0, -3.0);
        assert_ne!(a, b);
    }

    #[test]
    fn add_vector_and_point() {
        let a1 = Point::new(3.0, -2.0, 5.0);
        let a2 = Vector::new(-2.0, 3.0, 1.0);
        let result = a1 + a2;
        assert_eq!(result, Point::new(1.0, 1.0, 6.0));
    }

    #[test]
    fn point_sub() {
        let p1 = Point::new(3.0, 2.0, 1.0);
        let p2 = Point::new(5.0, 6.0, 7.0);
        let result = p1 - p2;
        assert_eq!(result, Vector::new(-2.0, -4.0, -6.0));
    }

    #[test]
    fn point_vector_sub() {
        let p1 = Point::new(3.0, 2.0, 1.0);
        let p2 = Vector::new(5.0, 6.0, 7.0);
        let result = p1 - p2;
        assert_eq!(result, Point::new(-2.0, -4.0, -6.0));
    }

    #[test]
    fn vector_sub() {
        let v1 = Vector::new(3.0, 2.0, 1.0);
        let v2 = Vector::new(5.0, 6.0, 7.0);
        let result = v1 - v2;
        assert_eq!(result, Vector::new(-2.0, -4.0, -6.0));
    }

    #[test]
    fn subtract_from_zero_vector_negates_vector() {
        let zero = Vector::new(0.0, 0.0, 0.0);
        let v = Vector::new(1.0, -2.0, 3.0);
        let result = zero - v;
        assert_eq!(result, Vector::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn negate_point_and_vector() {
        let a = Point::new(1.0, -2.0, 3.0);
        assert_eq!(-a, Point::new(-1.0, 2.0, -3.0));

        let b = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(-b, Vector::new(-1.0, 2.0, -3.0))
    }

    #[test]
    fn scalar_mul() {
        let a = Point::new(1.0, -2.0, 3.0);
        assert_eq!(a * 3.5, Point::new(3.5, -7.0, 10.5));

        let b = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(b * 3.5, Vector::new(3.5, -7.0, 10.5));
    }

    #[test]
    fn fraction_mul() {
        let a = Point::new(1.0, -2.0, 3.0);
        assert_eq!(a * 0.5, Point::new(0.5, -1.0, 1.5));

        let b = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(b * 0.5, Vector::new(0.5, -1.0, 1.5));
    }

    #[test]
    fn scalar_div() {
        let a = Point::new(1.0, -2.0, 3.0);
        assert_eq!(a / 2., Point::new(0.5, -1.0, 1.5));

        let b = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(b / 2., Vector::new(0.5, -1.0, 1.5));
    }

    #[test]
    fn unit_vector_magnitude() {
        let v = Vector::new(1.0, 0.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);

        let v = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);

        let v = Vector::new(0.0, 0.0, 1.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn non_unit_vector_magnitude() {
        let v = Vector::new(1.0, 2.0, 3.0);
        assert!(equal(v.magnitude(), 14.0_f64.sqrt()));

        let v = Vector::new(-1.0, -2.0, -3.0);
        assert!(equal(v.magnitude(), 14.0_f64.sqrt()));
    }

    #[test]
    fn normalize_vector() {
        let v = Vector::new(4.0, 0.0, 0.0);
        let norm = v.norm();
        assert_eq!(norm, Vector::new(1.0, 0.0, 0.0));

        let v = Vector::new(1.0, 2.0, 3.0);
        let norm = v.norm();
        assert!(equal(norm.x, 0.26726) && equal(norm.y, 0.53452) && equal(norm.z, 0.80178));
    }

    #[test]
    fn dot_product() {
        let a = Vector::new(1.0, 2.0, 3.0);
        let b = Vector::new(2.0, 3.0, 4.0);
        assert_eq!(a.dot(b), 20.0);
    }

    #[test]
    fn cross_product() {
        let a = Vector::new(1.0, 2.0, 3.0);
        let b = Vector::new(2.0, 3.0, 4.0);
        assert_eq!(a.cross(b), Vector::new(-1.0, 2.0, -1.0));
        assert_eq!(b.cross(a), Vector::new(1.0, -2.0, 1.0));
    }

    #[test]
    fn reflecting_vector_approaching_at_45_degrees() {
        let vector = Vector::new(1., -1., 0.);
        let normal = Vector::new(0., 1., 0.);
        let reflection = vector.reflect(normal);
        assert_eq!(reflection, Vector::new(1., 1., 0.));
    }

    #[test]
    fn reflecting_vector_off_slanted_surface() {
        let vector = Vector::new(0., -1., 0.);
        let val = 2.0_f64.sqrt() / 2.;
        let normal = Vector::new(val, val, 0.);
        let reflection = vector.reflect(normal);
        assert_eq!(reflection, Vector::new(1., 0., 0.));
    }
}
