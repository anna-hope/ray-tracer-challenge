use std::{fs::File, io::Write};

use ray_tracer::{canvas::Canvas, color::Color, tick, Environment, Projectile, Tuple};

fn main() {
    let mut projectile = Projectile {
        position: Tuple::point(0.0, 1.0, 0.0),
        velocity: Tuple::vector(1.0, 1.8, 0.0).norm() * 11.25,
    };
    let environment = Environment {
        gravity: Tuple::vector(0.0, -0.1, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };

    let color = Color::new(1.0, 0.0, 0.0);

    let mut canvas = Canvas::new(900, 550);

    loop {
        let x = projectile.position.x.round() as usize;
        let y = projectile.position.y.round() as usize;
        if let Some(y) = canvas.height.checked_sub(y) {
            canvas.write_pixel(x, y, color);
        }
        projectile = tick(&environment, &projectile);

        if projectile.position.y <= 0.0 {
            break;
        }
    }

    let ppm = canvas.to_ppm();
    let mut buffer = File::create("projectile.ppm").unwrap();
    buffer.write_all(ppm.as_bytes()).unwrap();
}
