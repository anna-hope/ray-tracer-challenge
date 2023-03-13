mod common;
mod tuple;

pub use tuple::Tuple;

pub struct Projectile {
    pub position: Tuple,
    pub velocity: Tuple,
}

pub struct Environment {
    pub gravity: Tuple,
    pub wind: Tuple,
}

pub fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    let position = proj.position + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    Projectile { position, velocity }
}
