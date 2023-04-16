use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use ray_tracer::prelude::*;

fn main() {
    let hex: Arc<dyn Shape> = Arc::new(group::Group::default());
    let hex_clone = Arc::clone(&hex);
    register_shape(hex_clone);

    for n in 0..=5 {
        let n = n as f64;
        let mut side: Arc<dyn Shape> = Arc::new(
            group::Group::default().with_transformation(Matrix::rotation_y(n * (PI / 3.))),
        );
        hex.add_child(&mut side);

        let mut corner: Arc<dyn Shape> = Arc::new(
            sphere::Sphere::default().with_transformation(
                Matrix::identity()
                    .scale(0.25, 0.25, 0.25)
                    .translate(0., 0., -1.),
            ),
        );

        let mut edge: Arc<dyn Shape> = Arc::new(cylinder::Cylinder::new(
            Matrix::identity()
                .scale(0.25, 1., 0.25)
                .rotate_z(-PI / 2.)
                .rotate_y(-PI / 6.)
                .translate(0., 0., -1.),
            Material::default(),
            0.,
            1.,
            false,
            None,
        ));

        side.add_child(&mut corner);
        side.add_child(&mut edge);
    }

    let light = Light::new(Tuple::point(0., 1., 1.), Color::white());
    let world = World::new(vec![hex], vec![light]);

    let from = Tuple::point(0., 4., 5.);
    let to = Tuple::point(0., 0., 0.);
    let up = Tuple::vector(0., 1., 0.);

    let camera_transformation = compute_view_transformation(from, to, up)
        .expect("Should be able to compute the camera transformation");
    let camera = Camera::new(500, 500, 0.5, camera_transformation);

    let canvas = camera
        .render(&world)
        .expect("Should be able to render the world");

    let output = canvas.to_ppm();
    let mut file = File::create("hexagon.ppm").expect("Should be able to create the output file");
    file.write_all(output.as_bytes())
        .expect("Should be able to write the output");

    println!("Wrote output to hexagon.ppm");
}
