use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::{common::equal, error::RayTracerError, Result};

pub enum TupleKind {
    Point,
    Vector,
}

#[derive(Debug, Clone, Copy)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Tuple {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    pub fn point(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 1.0 }
    }

    pub fn vector(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 0.0 }
    }

    pub fn kind(&self) -> TupleKind {
        if equal(self.w, 1.0) {
            TupleKind::Point
        } else {
            TupleKind::Vector
        }
    }

    pub fn magnitude(&self) -> f64 {
        let sum = self.x.powi(2) + self.y.powi(2) + self.z.powi(2) + self.w.powi(2);
        sum.sqrt()
    }

    pub fn norm(&self) -> Self {
        let magnitude = self.magnitude();
        Self {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
            w: self.w / magnitude,
        }
    }

    pub fn dot(&self, other: &Self) -> Result<f64> {
        // dotting non-vectors is a bug
        // the book says to use w (the last value) of a tuple
        // to tell if we've accidentally tried to do that
        // but thanks to Rust, we can do better --
        // return an actual error if that's what someone tried to do
        if !matches!(self.kind(), TupleKind::Vector) || !matches!(other.kind(), TupleKind::Vector) {
            return Err(RayTracerError::NonVectorTuple);
        }

        // the w (last component) for vectors is always 0, so no need to include it
        Ok(self.x * other.x + self.y * other.y + self.z * other.z)
    }

    pub fn cross(&self, other: &Self) -> Result<Self> {
        // ditto as with dot about returning error for non-vector operands
        if !matches!(self.kind(), TupleKind::Vector) || !matches!(other.kind(), TupleKind::Vector) {
            return Err(RayTracerError::NonVectorTuple);
        }

        // only implemented for the three-dimensional case,
        // since that's all we need (w is ignored since it's 0)
        Ok(Tuple::vector(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        ))
    }

    pub fn reflect(self, normal: Self) -> Result<Self> {
        if !matches!(self.kind(), TupleKind::Vector) || !matches!(normal.kind(), TupleKind::Vector)
        {
            return Err(RayTracerError::NonVectorTuple);
        }
        Ok(self - normal * 2. * self.dot(&normal)?)
    }
}

impl PartialEq for Tuple {
    fn eq(&self, other: &Self) -> bool {
        equal(self.x, other.x)
            && equal(self.y, other.y)
            && equal(self.z, other.z)
            && equal(self.w, other.w)
    }
}

impl Add for Tuple {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w,
        }
    }
}

impl Sub for Tuple {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w,
        }
    }
}

impl Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}

impl Mul<f64> for Tuple {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

impl Div<f64> for Tuple {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn point_is_point_and_not_vector() {
        let a = Tuple::new(4.3, -4.2, 3.1, 1.0);
        let kind = a.kind();
        assert!(matches!(kind, TupleKind::Point));
        assert!(!matches!(kind, TupleKind::Vector));
    }

    #[test]
    fn point_makes_point() {
        let a = Tuple::point(4.0, -4.0, 3.0);
        assert_eq!(a.w, 1.0);
    }

    #[test]
    fn vector_is_vector_and_not_point() {
        let a = Tuple::new(4.3, -4.2, 3.1, 0.0);
        let kind = a.kind();
        assert!(matches!(kind, TupleKind::Vector));
        assert!(!matches!(kind, TupleKind::Point));
    }

    #[test]
    fn vector_makes_vector() {
        let a = Tuple::vector(4.0, -4.0, 3.0);
        assert_eq!(a.w, 0.0);
    }

    #[test]
    fn two_tuples_equal() {
        let a = Tuple::point(1.0, 1.0, 1.0);
        let b = Tuple::point(1.0, 1.0, 1.0);
        assert_eq!(a, b);

        let a = Tuple::vector(1.0, 2.0, -3.0);
        let b = Tuple::vector(1.0, 2.0, -3.0);
        assert_eq!(a, b);
    }

    #[test]
    fn two_tuples_not_equal() {
        let a = Tuple::point(1.0, 1.0, -1.0);
        let b = Tuple::point(1.0, 1.0, 1.0);
        assert_ne!(a, b);

        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(1.0, 2.0, -3.0);
        assert_ne!(a, b);
    }

    #[test]
    fn add() {
        let a1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let a2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);
        let result = a1 + a2;
        assert_eq!(result, Tuple::new(1.0, 1.0, 6.0, 1.0));
    }

    #[test]
    fn point_sub() {
        let p1 = Tuple::point(3.0, 2.0, 1.0);
        let p2 = Tuple::point(5.0, 6.0, 7.0);
        let result = p1 - p2;
        assert_eq!(result, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn point_vector_sub() {
        let p1 = Tuple::point(3.0, 2.0, 1.0);
        let p2 = Tuple::vector(5.0, 6.0, 7.0);
        let result = p1 - p2;
        assert_eq!(result, Tuple::point(-2.0, -4.0, -6.0));
    }

    #[test]
    fn vector_sub() {
        let v1 = Tuple::vector(3.0, 2.0, 1.0);
        let v2 = Tuple::vector(5.0, 6.0, 7.0);
        let result = v1 - v2;
        assert_eq!(result, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn subtract_from_zero_vector_negates_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let v = Tuple::vector(1.0, -2.0, 3.0);
        let result = zero - v;
        assert_eq!(result, Tuple::vector(-1.0, 2.0, -3.0));
    }

    #[test]
    fn negate_tuple() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(-a, Tuple::new(-1.0, 2.0, -3.0, 4.0));
    }

    #[test]
    fn scalar_mul() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a * 3.5, Tuple::new(3.5, -7.0, 10.5, -14.0))
    }

    #[test]
    fn fraction_mul() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a * 0.5, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn scalar_div() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a / 2.0, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn unit_vector_magnitude() {
        let v = Tuple::vector(1.0, 0.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);

        let v = Tuple::vector(0.0, 1.0, 0.0);
        assert_eq!(v.magnitude(), 1.0);

        let v = Tuple::vector(0.0, 0.0, 1.0);
        assert_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn non_unit_vector_magnitude() {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        assert!(equal(v.magnitude(), 14.0_f64.sqrt()));

        let v = Tuple::vector(-1.0, -2.0, -3.0);
        assert!(equal(v.magnitude(), 14.0_f64.sqrt()));
    }

    #[test]
    fn normalize_vector() {
        let v = Tuple::vector(4.0, 0.0, 0.0);
        let norm = v.norm();
        assert_eq!(norm, Tuple::vector(1.0, 0.0, 0.0));

        let v = Tuple::vector(1.0, 2.0, 3.0);
        let norm = v.norm();
        assert!(equal(norm.x, 0.26726) && equal(norm.y, 0.53452) && equal(norm.z, 0.80178));
    }

    #[test]
    fn dot_product() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.dot(&b).unwrap(), 20.0);
    }

    #[test]
    fn non_vector_dot_gives_error() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::point(2.0, 3.0, 4.0);
        assert!(a.dot(&b).is_err());
    }

    #[test]
    fn cross_product() {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        assert_eq!(a.cross(&b).unwrap(), Tuple::vector(-1.0, 2.0, -1.0));
        assert_eq!(b.cross(&a).unwrap(), Tuple::vector(1.0, -2.0, 1.0));
    }

    #[test]
    fn reflecting_vector_approaching_at_45_degrees() {
        let vector = Tuple::vector(1., -1., 0.);
        let normal = Tuple::vector(0., 1., 0.);
        let reflection = vector.reflect(normal).unwrap();
        assert_eq!(reflection, Tuple::vector(1., 1., 0.));
    }

    #[test]
    fn reflecting_vector_off_slanted_surface() {
        let vector = Tuple::vector(0., -1., 0.);
        let val = 2.0_f64.sqrt() / 2.;
        let normal = Tuple::vector(val, val, 0.);
        let reflection = vector.reflect(normal).unwrap();
        assert_eq!(reflection, Tuple::vector(1., 0., 0.));
    }
}
