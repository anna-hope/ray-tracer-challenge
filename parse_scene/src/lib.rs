use serde_yaml::{Mapping, Sequence, Value};
use thiserror::Error;

use ray_tracer::{
    camera::Camera, light::Light, material::Material, shape::*,
    transformation::compute_view_transformation, Color, Matrix, Tuple,
};

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("YAML deserialization error: {0}")]
    Yaml(#[from] serde_yaml::Error),

    #[error("Ray tracer error: {0}")]
    RayTracer(#[from] ray_tracer::RayTracerError),

    #[error("Invalid YAML object: {0}")]
    InvalidYamlObject(String),

    #[error("Encountered unknown item: {0}")]
    UnknownItem(String),

    #[error("Invalid parameters for camera: {0}")]
    InvalidCameraParams(String),

    #[error("Invalid parameters for light: {0}")]
    InvalidLightParams(String),

    #[error("Invalid parameters for shape: {0}")]
    InvalidShapeParams(String),

    #[error("Invalid value")]
    InvalidValue,
}

pub type Result<T> = std::result::Result<T, ParseError>;

macro_rules! get_or_default_f64 {
    ($mapping:expr, $key:expr, $default:expr) => {
        $mapping
            .get($key)
            .and_then(|x| x.as_f64())
            .unwrap_or($default)
    };
}

#[derive(Debug)]
pub struct Scene {
    pub camera: Camera,
    pub lights: Vec<Light>,
    pub objects: Vec<Box<dyn Shape>>,
}

trait GetF64 {
    fn get_f64(&self) -> Result<f64>;
}

trait GetF64Vec {
    fn get_f64_vec(&self) -> Result<Vec<f64>>;
}

impl GetF64 for Value {
    fn get_f64(&self) -> Result<f64> {
        self.as_f64().ok_or(ParseError::InvalidValue)
    }
}

impl GetF64Vec for Value {
    fn get_f64_vec(&self) -> Result<Vec<f64>> {
        self.as_sequence()
            .ok_or(ParseError::InvalidValue)?
            .iter()
            .map(|f| f.get_f64())
            .collect::<Result<Vec<f64>>>()
    }
}

fn construct_camera(description: &Mapping) -> Result<Camera> {
    let width = description
        .get("width")
        .ok_or(ParseError::InvalidCameraParams(
            "Missing camera width".to_string(),
        ))?
        .as_u64()
        .ok_or(ParseError::InvalidCameraParams(
            "Camera width must be an integer".to_string(),
        ))? as usize;

    let height = description
        .get("height")
        .ok_or(ParseError::InvalidCameraParams(
            "Missing camera height".to_string(),
        ))?
        .as_u64()
        .ok_or(ParseError::InvalidCameraParams(
            "Camera height must be an integer".to_string(),
        ))? as usize;

    let field_of_view = description
        .get("field-of-view")
        .ok_or(ParseError::InvalidCameraParams(
            "Missing field-of-view".to_string(),
        ))?
        .get_f64()?;

    let from_values = description
        .get("from")
        .ok_or(ParseError::InvalidCameraParams(
            "Missing 'from'".to_string(),
        ))?
        .get_f64_vec()?;

    let to_values = description
        .get("to")
        .ok_or(ParseError::InvalidCameraParams("Missing 'to'".to_string()))?
        .get_f64_vec()?;

    let up_values = description
        .get("up")
        .ok_or(ParseError::InvalidCameraParams("Missing 'up'".to_string()))?
        .get_f64_vec()?;

    let from = Tuple::point(from_values[0], from_values[1], from_values[2]);
    let to = Tuple::point(to_values[0], to_values[1], to_values[2]);
    let up = Tuple::vector(up_values[0], up_values[1], up_values[2]);

    let transformation = compute_view_transformation(from, to, up)?;
    let camera = Camera::new(width, height, field_of_view).with_transformation(transformation);

    Ok(camera)
}

fn construct_light(description: &Mapping) -> Result<Light> {
    let position_values = description
        .get("at")
        .ok_or(ParseError::InvalidLightParams("Missing 'at'".to_string()))?
        .get_f64_vec()?;

    let intensity_values = description
        .get("intensity")
        .ok_or(ParseError::InvalidLightParams(
            "Missing 'intensity'".to_string(),
        ))?
        .get_f64_vec()?;

    let position = Tuple::point(position_values[0], position_values[1], position_values[2]);
    let intensity = Color::new(
        intensity_values[0],
        intensity_values[1],
        intensity_values[2],
    );

    Ok(Light::new(position, intensity))
}

fn parse_material(description: &Value) -> Material {
    let color_values = description
        .get("color")
        .and_then(|x| x.get_f64_vec().ok())
        .unwrap_or_else(|| vec![1., 1., 1.]);

    let color = Color::new(color_values[0], color_values[1], color_values[2]);

    let ambient = get_or_default_f64!(description, "ambient", 0.1);
    let diffuse = get_or_default_f64!(description, "diffuse", 0.9);
    let specular = get_or_default_f64!(description, "specular", 0.9);
    let shininess = get_or_default_f64!(description, "shininess", 200.);
    let reflectivity = get_or_default_f64!(description, "reflectivity", 0.);
    let transparency = get_or_default_f64!(description, "transparency", 0.);
    let refractive_index = get_or_default_f64!(description, "refractive-index", 1.);
    let casts_shadow = description
        .get("casts-shadow")
        .and_then(|x| x.as_bool())
        .unwrap_or(true);

    Material {
        ambient,
        diffuse,
        color,
        specular,
        shininess,
        reflectivity,
        transparency,
        refractive_index,
        casts_shadow,
        ..Default::default()
    }
}

pub fn construct_object(object_type: &str, description: &Mapping) -> Result<Box<dyn Shape>> {
    let material = if let Some(material_description) = description.get("material") {
        parse_material(material_description)
    } else {
        Material::default()
    };

    let mut transformation = Matrix::identity();

    if let Some(transformations) = description.get("transform").and_then(|x| x.as_sequence()) {
        for value in transformations {
            let key = value
                .get(0)
                .ok_or(ParseError::InvalidValue)?
                .as_str()
                .ok_or(ParseError::InvalidValue)?;

            let value = value.get(1).ok_or(ParseError::InvalidValue)?;

            if value.is_sequence() {
                let value = value.get_f64_vec()?;

                match key {
                    "scale" => {
                        transformation = transformation.scale(value[0], value[1], value[2]);
                    }
                    "translate" => {
                        transformation = transformation.translate(value[0], value[1], value[2]);
                    }
                    &_ => unimplemented!(),
                }
            } else {
                let radians = value.get_f64()?;

                match key {
                    "rotate-x" => {
                        transformation = transformation.rotate_x(radians);
                    }
                    "rotate-y" => {
                        transformation = transformation.rotate_y(radians);
                    }
                    "rotate-z" => {
                        transformation = transformation.rotate_z(radians);
                    }
                    _ => unimplemented!(),
                }
            }
        }
    }

    let object: Box<dyn Shape> = match object_type {
        "sphere" => Box::new(sphere::Sphere::new(transformation, material)),
        "plane" => Box::new(plane::Plane::new(transformation, material)),
        "cube" => Box::new(cube::Cube::new(transformation, material)),
        _ => unimplemented!(),
    };

    Ok(object)
}

pub fn parse_scene(input: &str) -> Result<Scene> {
    let sequence = serde_yaml::from_str::<Sequence>(input)?;

    let mut camera: Option<Camera> = None;
    let mut lights = vec![];
    let mut objects: Vec<Box<dyn Shape>> = vec![];

    for item in sequence {
        if let Some(mapping) = item.as_mapping() {
            if let Some(Value::String(item_type)) = mapping.get("add") {
                match item_type.as_str() {
                    "camera" => {
                        camera = Some(construct_camera(mapping)?);
                    }
                    "light" => {
                        let light = construct_light(mapping)?;
                        lights.push(light);
                    }
                    "sphere" | "plane" | "cube" => {
                        let object = construct_object(item_type.as_str(), mapping)?;
                        objects.push(object);
                    }
                    _ => return Err(ParseError::UnknownItem(item_type.to_owned())),
                }
            }
        } else {
            return Err(ParseError::InvalidYamlObject(format!(
                "Expected mapping, found {item:#?}"
            )));
        }
    }

    let camera = camera.expect("The scene must define a camera!");

    let scene = Scene {
        camera,
        lights,
        objects,
    };
    Ok(scene)
}
