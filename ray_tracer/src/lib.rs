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
pub use matrix::Matrix;
pub use tuple::Tuple;

#[derive(Debug, PartialEq)]
pub enum SceneObjectType {
    Sphere,
}

pub trait SceneObject: intersection::Intersect {
    /// Compute the normal vector at the world point.
    fn normal_at(&self, world_point: Tuple) -> Tuple;

    /// Get object id.
    fn id(&self) -> usize;

    /// Get object type.
    fn object_type(&self) -> SceneObjectType;

    fn material(&self) -> material::Material;
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
