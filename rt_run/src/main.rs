use ray_tracer::{tick, Environment, Projectile, Tuple};

fn main() {
    let mut projectile = Projectile {
        position: Tuple::point(0.0, 1.0, 0.0),
        velocity: Tuple::vector(1.0, 1.0, 0.0).norm(),
    };
    let environment = Environment {
        gravity: Tuple::vector(0.0, -0.1, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };

    let mut num_ticks: u32 = 0;
    loop {
        projectile = tick(&environment, &projectile);
        println!("Projectile position: {:?}", projectile.position);
        num_ticks += 1;

        if projectile.position.y <= 0.0 {
            break;
        }
    }
    println!("Took {num_ticks} ticks for projectile to hit the ground.");
}
