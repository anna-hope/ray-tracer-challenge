pub mod material;
pub mod canvas;
pub mod color;
mod common;
pub mod intersection;
pub mod light;
mod matrix;
pub mod sphere;
mod tuple;

use std::fmt::Debug;

pub use matrix::Matrix;
pub use tuple::Tuple;

#[derive(Debug, PartialEq)]
pub enum SceneObjectType {
    Sphere,
}

pub trait Id {
    fn id(&self) -> usize;
}

pub trait ObjectType {
    fn object_type(&self) -> SceneObjectType;
}

pub trait SceneObject: intersection::Intersect + ObjectType + Id {}

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
