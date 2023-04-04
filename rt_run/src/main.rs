use std::{f64::consts::PI, fs::File, io::Write};

use ray_tracer::{
    camera,
    color::Color,
    light::PointLight,
    material::Material,
    pattern::*,
    shape::{plane::Plane, sphere::Sphere},
    transformation,
    world::World,
    Matrix, Tuple,
};

fn main() {
    let solid_pattern1 = Box::new(solid::SolidPattern::new(Color::new(0.8, 0.8, 1.0)));
    let solid_pattern2 = Box::new(solid::SolidPattern::new(Color::from_rgb(32, 178, 170)));

    let blended_pattern1 =
        blended::BlendedPattern::new(solid_pattern1, solid_pattern2, Matrix::identity());

    let solid_pattern1 = Box::new(solid::SolidPattern::new(Color::from_rgb(100, 40, 77)));
    let solid_pattern2 = Box::new(solid::SolidPattern::new(Color::from_rgb(7, 200, 90)));

    let blended_pattern2 =
        blended::BlendedPattern::new(solid_pattern1, solid_pattern2, Matrix::identity());

    let floor_pattern1 = checker::CheckerPattern::new(
        Box::new(blended_pattern1),
        Box::new(blended_pattern2),
        Matrix::identity(),
    );

    let solid_pattern1 = solid::SolidPattern::new(Color::from_rgb(100, 200, 50));
    let solid_pattern2 = solid::SolidPattern::new(Color::from_rgb(13, 5, 70));

    let floor_pattern2 = checker::CheckerPattern::new(
        Box::new(solid_pattern1),
        Box::new(solid_pattern2),
        Matrix::identity(),
    );

    let floor_pattern = checker::CheckerPattern {
        a: Box::new(floor_pattern1),
        b: Box::new(floor_pattern2),
        ..Default::default()
    };
    let floor_material = Material {
        pattern: Some(Box::new(floor_pattern)),
        reflectivity: 0.6,
        ..Default::default()
    };
    let floor = Plane::new()
        .with_material(floor_material)
        .with_transformation(Matrix::scaling(3., 3., 3.));

    let wall_transformation = Matrix::identity().rotate_x(PI / 2.).translate(0., 0., 3.);
    let wall_material = Material {
        reflectivity: 1.,
        ..Default::default()
    };
    let wall = Plane::new()
        .with_transformation(wall_transformation)
        .with_material(wall_material);

    let middle_sphere_transformation = Matrix::translation(-0.5, 1., 0.5);
    // let middle_sphere_pattern = Box::new(SolidPattern::new(Color::from_rgb(75, 0, 130)));

    let middle_sphere_material = Material {
        color: Color::from_rgb(75, 0, 130),
        transparency: 0.9,
        reflectivity: 1.,
        diffuse: 0.1,
        specular: 1.,
        shininess: 300.,
        ambient: 0.1,
        pattern: None,
        casts_shadow: false,
        ..Default::default()
    };
    let middle_sphere = Sphere::default()
        .with_material(middle_sphere_material)
        .with_transformation(middle_sphere_transformation);

    let right_sphere_transformation = Matrix::identity()
        .scale(0.5, 0.5, 0.5)
        .translate(1.5, 0.5, -0.5);
    let right_sphere_material = Material {
        color: Color::new(0.5, 1., 0.1),
        diffuse: 0.7,
        specular: 0.3,
        reflectivity: 0.5,
        pattern: Some(Box::new(gradient::GradientPattern {
            a: Color::from_rgb(218, 112, 214),
            b: Color::from_rgb(75, 0, 130),
            transformation: Matrix::scaling(0.5, 0.5, 0.5),
        })),
        ..Default::default()
    };
    let right_sphere = Sphere::default()
        .with_transformation(right_sphere_transformation)
        .with_material(right_sphere_material);

    let left_sphere_subpattern = radial_gradient::RadialGradientPattern::new(
        Color::new(0.1, 0., 0.9),
        Color::new(1., 0., 1.),
        Matrix::identity(),
    );

    let left_sphere_pattern = perturbed::PerturbedPattern::new(
        Box::new(left_sphere_subpattern),
        Matrix::scaling(0.2, 0.2, 0.2),
    );

    let left_sphere_transformation = Matrix::scaling(0.33, 0.33, 0.33).translate(-1.5, 0.33, -0.75);
    let left_sphere_material = Material {
        color: Color::new(1., 0.8, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::new(left_sphere_pattern)),
        ..Default::default()
    };
    let left_sphere = Sphere::default()
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
        Tuple::point(0., 5., -7.),
        Tuple::point(0., 2., 0.),
        Tuple::vector(0., 1., 0.),
    )
    .expect("Should compute camera transformation");
    let camera = camera::Camera::new(1000, 500, PI / 3.).with_transformation(camera_transformation);
    let canvas = camera.render(&world).expect("Should render the image");

    let ppm = canvas.to_ppm();
    let mut buffer = File::create("world.ppm").unwrap();
    buffer.write_all(ppm.as_bytes()).unwrap();
}
