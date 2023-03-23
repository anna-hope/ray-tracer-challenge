pub mod camera;
pub mod canvas;
pub mod color;
mod error;
pub mod intersection;
pub mod light;
pub mod material;
mod matrix;
pub mod sphere;
pub mod transformation;
mod tuple;
pub mod world;

use std::fmt::Debug;

pub use color::Color;
use error::RayTracerError;
pub use intersection::Ray;
pub use matrix::Matrix;
pub use tuple::Tuple;

type Result<T> = std::result::Result<T, RayTracerError>;

const EPSILON: f64 = 1e-5;

#[derive(Debug, PartialEq)]
pub enum ShapeType {
    Sphere,
}

pub trait Shape: intersection::Intersect {
    /// Computes the normal vector at the world point.
    fn normal_at(&self, world_point: Tuple) -> Tuple;

    /// Gets object id.
    fn id(&self) -> usize;

    /// Gets object type.
    fn shape_type(&self) -> ShapeType;

    fn material(&self) -> material::Material;

    /// Gets an intersection with an arbitrary t for this object.
    /// This is needed primarily for testing, because we can't construct
    /// an intersection with a boxed trait object due to type incompatibility.
    fn arbitrary_intersection(&self, t: f64) -> intersection::Intersection;

    /// Sets the object material to the given material.
    /// Needed primarily for testing.
    fn set_material(&mut self, material: material::Material);
}

impl PartialEq for dyn Shape {
    fn eq(&self, other: &Self) -> bool {
        self.shape_type() == other.shape_type() && self.id() == other.id()
    }
}

impl Debug for dyn Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SceneObject")
            .field("type", &self.shape_type())
            .field("id", &self.id())
            .finish()
    }
}

pub fn equal(a: f64, b: f64) -> bool {
    let c = a - b;
    c.abs() < EPSILON
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn two_floats_equal() {
        let a = 1.0;
        let b = 1.00000001;
        assert!(equal(a, b));
    }
}
