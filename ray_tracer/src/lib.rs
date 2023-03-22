pub mod canvas;
pub mod color;
mod common;
pub mod intersection;
pub mod light;
pub mod material;
mod matrix;
pub mod sphere;
mod tuple;
pub mod world;

use std::fmt::Debug;

pub use color::Color;
pub use intersection::Ray;
pub use matrix::Matrix;
pub use tuple::Tuple;

#[derive(Debug, PartialEq)]
pub enum SceneObjectType {
    Sphere,
}

pub trait SceneObject: intersection::Intersect {
    /// Computes the normal vector at the world point.
    fn normal_at(&self, world_point: Tuple) -> Tuple;

    /// Gets object id.
    fn id(&self) -> usize;

    /// Gets object type.
    fn object_type(&self) -> SceneObjectType;

    fn material(&self) -> material::Material;

    /// Gets an intersection with an arbitrary t for this object.
    /// This is needed primarily for testing, because we can't construct
    /// an intersection with a boxed trait object due to type incompatibility.
    fn arbitrary_intersection(&self, t: f64) -> intersection::Intersection;

    /// Sets the object material to the given material.
    /// Needed primarily for testing.
    fn set_material(&mut self, material: material::Material);
}

impl PartialEq for dyn SceneObject {
    fn eq(&self, other: &Self) -> bool {
        self.object_type() == other.object_type() && self.id() == other.id()
    }
}

impl Debug for dyn SceneObject {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SceneObject")
            .field("type", &self.object_type())
            .field("id", &self.id())
            .finish()
    }
}
