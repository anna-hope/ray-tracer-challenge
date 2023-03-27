use std::{f64::consts::PI, fs::File, io::Write};

use ray_tracer::{
    camera,
    color::Color,
    light::PointLight,
    material::Material,
    pattern,
    shape::{plane::Plane, sphere::Sphere},
    transformation,
    world::World,
    Matrix, Tuple,
};

fn main() {
    let floor_pattern = pattern::CheckerPattern {
        a: Color::new(0.5, 0.2, 0.8),
        b: Color::new(0.9, 0.5, 0.1),
        ..Default::default() // transformation: Matrix::scaling(3., 3., 3.),
    };
    let floor_material = Material {
        pattern: Some(Box::new(floor_pattern.clone())),
        ..Default::default()
    };
    let floor = Plane::new().with_material(floor_material);
    // .with_transformation(Matrix::scaling(3., 3., 3.));

    let wall_transformation = Matrix::identity().rotate_x(PI / 2.).translate(0., 0., 3.);
    let wall = Plane::new().with_transformation(wall_transformation);

    let middle_sphere_transformation = Matrix::translation(-0.5, 1., 0.5);
    let middle_sphere_pattern = pattern::CheckerPattern {
        a: Color::new(0.1, 0., 0.9),
        b: Color::white(),
        ..Default::default()
    }
    .with_transformation(Matrix::scaling(0.5, 1.0, 0.5));

    let middle_sphere_material = Material {
        color: Color::new(0.1, 1., 0.5),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::new(middle_sphere_pattern)),
        ..Default::default()
    };
    let middle_sphere = Sphere::new()
        .with_material(middle_sphere_material)
        .with_transformation(middle_sphere_transformation);

    let right_sphere_transformation = Matrix::identity()
        .scale(0.5, 0.5, 0.5)
        .translate(1.5, 0.5, -0.5);
    let right_sphere_material = Material {
        color: Color::new(0.5, 1., 0.1),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::new(pattern::StripePattern {
            a: Color::new(0.1, 0., 0.9),
            b: Color::white(),
            ..Default::default()
        })),
        ..Default::default()
    };
    let right_sphere = Sphere::new()
        .with_transformation(right_sphere_transformation)
        .with_material(right_sphere_material);

    let left_sphere_transformation = Matrix::identity()
        .scale(0.33, 0.33, 0.33)
        .translate(-1.5, 0.33, -0.75);
    let left_sphere_material = Material {
        color: Color::new(1., 0.8, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::new(pattern::StripePattern {
            a: Color::new(0.1, 0., 0.9),
            b: Color::white(),
            ..Default::default()
        })),
        ..Default::default()
    };
    let left_sphere = Sphere::new()
        .with_transformation(left_sphere_transformation)
        .with_material(left_sphere_material);

    let world = World {
        light: Some(PointLight::new(
            Tuple::point(-10., 10., -10.),
            Color::white(),
        )),
        objects: vec![
            Box::new(floor),
            Box::new(wall),
            Box::new(left_sphere),
            Box::new(middle_sphere),
            Box::new(right_sphere),
        ],
    };

    let camera_transformation = transformation::compute_view_transformation(
        Tuple::point(0., 1.5, -5.),
        Tuple::point(0., 1., 0.),
        Tuple::vector(0., 1., 0.),
    )
    .expect("Should compute camera transformation");
    let camera = camera::Camera::new(1000, 500, PI / 3.).with_transformation(camera_transformation);
    let canvas = camera.render(&world).expect("Should render the image");

    let ppm = canvas.to_ppm();
    let mut buffer = File::create("world.ppm").unwrap();
    buffer.write_all(ppm.as_bytes()).unwrap();
}
