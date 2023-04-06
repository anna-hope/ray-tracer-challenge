pub mod camera;
pub mod canvas;
pub mod color;
mod error;
pub mod intersection;
pub mod light;
pub mod material;
mod matrix;
pub mod pattern;
pub mod shape;
pub mod transformation;
mod tuple;
pub mod world;

pub use color::Color;
pub use error::RayTracerError;
pub use intersection::Ray;
pub use matrix::Matrix;
pub use tuple::Tuple;

type Result<T> = std::result::Result<T, RayTracerError>;

const EPSILON: f64 = 1e-5;

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
