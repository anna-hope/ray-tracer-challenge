use thiserror::Error;

#[derive(Debug, Error)]
pub enum RayTracerError {
    #[error("Attempt to use a non-vector tuple in a vector-only context")]
    NonVectorTuple,
}
