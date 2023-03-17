#![allow(dead_code)]

use std::{f64::consts::PI, fs::File, io::Write};

use ray_tracer::{canvas::Canvas, color::Color, Matrix, Tuple};

struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    let position = proj.position + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    Projectile { position, velocity }
}

fn main() {
    let color = Color::new(1.0, 1.0, 1.0);

    let canvas_width = 300;
    let mut canvas = Canvas::new(canvas_width, canvas_width);
    let center_x = canvas_width / 2;
    let center_y = center_x;
    let clock_radius = (canvas_width as f64) * 3. / 8.;

    let twelve = Tuple::point(0., 0., 1.);

    for num in 0..12 {
        let rotation = Matrix::rotation_y((num as f64) * PI / 6.);
        let point = rotation * twelve;

        // need to add center_x and center_y *before* we convert f64 to usize
        // because when a negative f64 is converted to usize, it turns to 0
        let canvas_x = (point.x * clock_radius + center_x as f64).round() as usize;
        let canvas_y = (point.z * clock_radius + center_y as f64).round() as usize;

        canvas.write_pixel(canvas_x, canvas_y, color);
    }

    let ppm = canvas.to_ppm();
    let mut buffer = File::create("clock.ppm").unwrap();
    buffer.write_all(ppm.as_bytes()).unwrap();
}
