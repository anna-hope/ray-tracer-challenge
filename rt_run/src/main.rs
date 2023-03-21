#![allow(dead_code)]

use std::{f64::consts::PI, fs::File, io::Write};

use indicatif::{ProgressBar, ProgressStyle};

use ray_tracer::{
    canvas::Canvas,
    color::Color,
    intersection::{hit, Intersect, Ray},
    sphere::Sphere,
    Matrix, Tuple,
};

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
    println!("Rendering...");
    let ray_origin = Tuple::point(0., 0., -5.);
    let wall_z = 10.;

    let canvas_width = 500;
    let wall_size = 7.;
    let pixel_size = wall_size / canvas_width as f64;
    let half = wall_size / 2.;

    let mut canvas = Canvas::new(canvas_width, canvas_width);
    let color = Color::new(1., 0., 0.);

    let transformation = Matrix::identity()
        .shear(1., 0., 0., 0., 0., 0.)
        .scale(0.5, 1., 1.);
    let shape = Sphere::new().with_transformation(transformation);

    let progress_bar = ProgressBar::new(canvas_width as u64);
    let progress_style =
        ProgressStyle::with_template("[{elapsed}] {eta} {bar} {pos:>7}/{len:7} {percent}% {msg}")
            .unwrap();
    progress_bar.set_style(progress_style);

    // for each row of pixels in the canvas
    for y in 0..canvas_width - 1 {
        progress_bar.inc(1);
        // compute the world y coordinates
        let world_y = half - (pixel_size * y as f64);

        for x in 0..canvas_width - 1 {
            let world_x = -half + (pixel_size * x as f64);

            // describe the point on the wall that the ray will target
            let position = Tuple::point(world_x, world_y, wall_z);
            let normalized_direction = (position - ray_origin).norm();
            let ray = Ray::new(ray_origin, normalized_direction);

            let xs = shape.intersect(&ray);
            if hit(&xs).is_some() {
                canvas.write_pixel(x, y, color);
            }
        }
    }

    let ppm = canvas.to_ppm();
    let mut buffer = File::create("sphere.ppm").unwrap();
    buffer.write_all(ppm.as_bytes()).unwrap();
}
