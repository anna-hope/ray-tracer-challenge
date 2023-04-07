use std::{
    env,
    fs::File,
    io::{Read, Write},
};

use parse_scene::parse_scene;
use ray_tracer::world::World;

fn main() {
    let mut args = env::args();

    let yaml_filename = args
        .nth(1)
        .expect("Must provide a path to a YAML file with a scene description.");
    let mut file = File::open(yaml_filename).expect("The file should exist and be openable");
    let mut data = String::new();
    file.read_to_string(&mut data)
        .expect("The file should be readable to a string");

    let scene = parse_scene(data.as_str()).expect("Should be able to parse scene YAML");

    let world = World::new(scene.objects, scene.lights.get(0).cloned());

    let canvas = scene
        .camera
        .render(&world)
        .expect("Should be able to render the scene");

    let ppm = canvas.to_ppm();
    let mut buffer = File::create("scene.ppm").unwrap();
    buffer.write_all(ppm.as_bytes()).unwrap();
}
