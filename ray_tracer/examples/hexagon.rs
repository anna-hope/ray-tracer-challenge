use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use ray_tracer::prelude::*;

fn main() {
    let hex = Arc::new(group::Group::default());
    let hex_clone = Arc::clone(&hex) as Arc<dyn Shape>;
    register_shape(hex_clone);

    for n in 0..=5 {
        let n = n as f64;
        let side: Arc<dyn Shape> = Arc::new(
            group::Group::default().with_transformation(Matrix::rotation_y(n * (PI / 3.))),
        );
        hex.add_child(&side);

        let corner: Arc<dyn Shape> = Arc::new(
            sphere::Sphere::default().with_transformation(
                Matrix::identity()
                    .scale(0.25, 0.25, 0.25)
                    .translate(0., 0., -1.),
            ),
        );

        let edge: Arc<dyn Shape> = Arc::new(cylinder::Cylinder::new(
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

        let side = side.as_group().unwrap();
        side.add_child(&corner);
        side.add_child(&edge);
    }

    let light = PointLight::new(Point::new(0., 1., 1.), Color::white());
    let world = World::new(vec![hex], vec![light]);

    let from = Point::new(0., 4., 5.);
    let to = Point::new(0., 0., 0.);
    let up = Vector::new(0., 1., 0.);

    let camera_transformation = compute_view_transformation(from, to, up);
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
