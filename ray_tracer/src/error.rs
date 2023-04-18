use thiserror::Error;

#[derive(Debug, Error)]
pub enum RayTracerError {
    #[error("Attempted to calculate an inverse of a non-invertible matrix")]
    NonInvertibleMatrix,
}
