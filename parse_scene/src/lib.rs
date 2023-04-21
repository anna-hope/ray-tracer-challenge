use std::sync::Arc;

use serde_yaml::{Mapping, Sequence, Value};
use thiserror::Error;

use ray_tracer::prelude::*;

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

    #[error("The scene is missing camera defintion")]
    MissingCamera,

    #[error("Invalid parameters for camera: {0}")]
    InvalidCameraParams(String),

    #[error("Invalid parameters for light: {0}")]
    InvalidLightParams(String),

    #[error("Invalid parameters for shape: {0}")]
    InvalidShapeParams(String),

    #[error("Invalid parameters for a constant definition: {0}")]
    InvalidDefinitionParams(String),

    #[error("Invalid parameters for transformation: {0}")]
    InvalidTransformationParams(String),

    #[error("Invalid value: {0}")]
    InvalidValue(String),

    #[error("Invalid parameters for group: {0}")]
    InvalidGroupParams(String),
}

type Result<T> = std::result::Result<T, ParseError>;

macro_rules! get_or_default_f64 {
    ($mapping:expr, $key:expr, $default:expr) => {
        $mapping
            .get($key)
            .and_then(|x| x.as_f64())
            .unwrap_or($default)
    };
}

#[derive(Debug, Clone, Copy)]
enum Transformation {
    Scale(f64, f64, f64),
    Translate(f64, f64, f64),
    RotateX(f64),
    RotateY(f64),
    RotateZ(f64),
}

#[derive(Debug)]
enum DefinitionValue {
    Transformations(Vec<Transformation>),
    Material(Material),
}

#[derive(Debug)]
struct Definition {
    name: String,
    value: DefinitionValue,
}

#[derive(Debug)]
pub struct Scene {
    pub camera: Camera,
    pub lights: Vec<Light>,
    pub objects: Vec<Arc<dyn Shape>>,
}

trait GetF64 {
    fn get_f64(&self) -> Result<f64>;
}

trait GetF64Vec {
    fn get_f64_vec(&self) -> Result<Vec<f64>>;
}

impl GetF64 for Value {
    fn get_f64(&self) -> Result<f64> {
        self.as_f64()
            .ok_or(ParseError::InvalidValue(format!("{self:#?}")))
    }
}

impl GetF64Vec for Value {
    fn get_f64_vec(&self) -> Result<Vec<f64>> {
        self.as_sequence()
            .ok_or(ParseError::InvalidValue(format!("{self:#?}")))?
            .iter()
            .map(|f| f.get_f64())
            .collect::<Result<Vec<f64>>>()
    }
}

fn parse_transformation(item: &Value) -> Result<Transformation> {
    // each transformation is described as an array of 2+ elements
    // where the first element is the string key (name) of the transformation
    // and the following elements are its numerical parameters

    let transformation_description = item
        .as_sequence()
        .ok_or(ParseError::InvalidValue(format!("{item:#?}")))?;
    let transformation_name = transformation_description
        .get(0)
        .ok_or(ParseError::InvalidTransformationParams(
            "Missing transformation name".to_string(),
        ))?
        .as_str()
        .ok_or(ParseError::InvalidValue(format!("{item:#?}")))?;

    // transformation names can be 'scale', 'translate', 'rotate-x', 'rotate-y', 'rotate-z'
    // 'scale' and 'translate' take 3 parameters and are described as [name, x, y, z]
    // 'rotate-x', 'rotate-y', 'rotate-z' take 1 parameter and are described as [name, angle]

    match transformation_name {
        "scale" => {
            let params = transformation_description
                .get(1..)
                .ok_or(ParseError::InvalidTransformationParams(
                    "Missing transformation parameters".to_string(),
                ))?
                .iter()
                .map(|f| f.get_f64())
                .collect::<Result<Vec<f64>>>()?;
            Ok(Transformation::Scale(params[0], params[1], params[2]))
        }
        "translate" => {
            let params = transformation_description
                .get(1..)
                .ok_or(ParseError::InvalidTransformationParams(
                    "Missing transformation parameters".to_string(),
                ))?
                .iter()
                .map(|f| f.get_f64())
                .collect::<Result<Vec<f64>>>()?;
            Ok(Transformation::Translate(params[0], params[1], params[2]))
        }
        "rotate-x" => {
            let angle = transformation_description
                .get(1)
                .ok_or(ParseError::InvalidTransformationParams(
                    "Missing transformation parameters".to_string(),
                ))?
                .get_f64()?;
            Ok(Transformation::RotateX(angle))
        }
        "rotate-y" => {
            let angle = transformation_description
                .get(1)
                .ok_or(ParseError::InvalidTransformationParams(
                    "Missing transformation parameters".to_string(),
                ))?
                .get_f64()?;
            Ok(Transformation::RotateY(angle))
        }
        "rotate-z" => {
            let angle = transformation_description
                .get(1)
                .ok_or(ParseError::InvalidTransformationParams(
                    "Missing transformation parameters".to_string(),
                ))?
                .get_f64()?;
            Ok(Transformation::RotateZ(angle))
        }
        _ => Err(ParseError::InvalidTransformationParams(format!(
            "Unknown transformation name: {transformation_name}"
        ))),
    }
}

fn parse_definition(description: &Mapping, definitions: &[Definition]) -> Result<Definition> {
    let name = description
        .get("define")
        .ok_or(ParseError::InvalidDefinitionParams(
            "Missing 'define' field".to_string(),
        ))?
        .as_str()
        .ok_or(ParseError::InvalidValue(format!("{description:#?}")))?
        .to_string();

    // a definition can either define a material
    // or a series of transformations
    // in both cases, it can extend another definition
    // if it extends a material definition, it can override its fields
    // if it extends a transformation definition, it can append to its list of transformations

    // this definition is a material definition if its 'value' field is a mapping

    let definition_value = description
        .get("value")
        .ok_or(ParseError::InvalidDefinitionParams(
            "Missing 'value' field".to_string(),
        ))?;

    if definition_value.is_mapping() {
        let material = if let Some(extend_value) = description.get("extend") {
            let value_str = extend_value
                .as_str()
                .ok_or(ParseError::InvalidValue(format!("{extend_value:#?}")))?;

            // get the original definition that this one is extending
            let original_definition = definitions.iter().find(|d| d.name == value_str).ok_or(
                ParseError::InvalidDefinitionParams(format!("Definition {} not found", value_str)),
            )?;
            let original_material = match &original_definition.value {
                DefinitionValue::Material(m) => m,
                _ => {
                    return Err(ParseError::InvalidDefinitionParams(
                        "Cannot extend a definition that does not have a material".to_string(),
                    ))
                }
            };

            // clone the material so that we can override the fields
            // with the values specified in the definition
            let mut material = original_material.clone();

            let new_material = parse_material(definition_value);

            // override the field of material with what is specified in the values
            if definition_value.get("color").is_some() {
                material.color = new_material.color;
            }

            if definition_value.get("ambient").is_some() {
                material.ambient = new_material.ambient;
            }

            if definition_value.get("diffuse").is_some() {
                material.diffuse = new_material.diffuse;
            }

            if definition_value.get("specular").is_some() {
                material.specular = new_material.specular;
            }

            if definition_value.get("shininess").is_some() {
                material.shininess = new_material.shininess;
            }

            if definition_value.get("reflective").is_some() {
                material.reflectivity = new_material.reflectivity;
            }

            if definition_value.get("refractive-index").is_some() {
                material.refractive_index = new_material.refractive_index;
            }

            if definition_value.get("transparency").is_some() {
                material.transparency = new_material.transparency;
            }

            if definition_value.get("casts-shadow").is_some() {
                material.casts_shadow = new_material.casts_shadow;
            }

            material
        } else {
            let material_spec =
                description
                    .get("value")
                    .ok_or(ParseError::InvalidDefinitionParams(
                        "Missing 'value'".to_string(),
                    ))?;
            parse_material(material_spec)
        };

        Ok(Definition {
            name,
            value: DefinitionValue::Material(material),
        })
    } else if let Some(seq) = definition_value.as_sequence() {
        // this extends other transformations if any of the elements are strings
        let mut transformations = Vec::new();

        for item in seq {
            if let Some(original_definition) = item.as_str() {
                // find the original definition
                let original_definition = definitions
                    .iter()
                    .find(|d| d.name == original_definition)
                    .ok_or(ParseError::InvalidDefinitionParams(format!(
                        "Definition {} not found",
                        original_definition
                    )))?;

                // get its transformations
                let mut original_transformations = match &original_definition.value {
                    DefinitionValue::Transformations(ts) => ts.to_owned(),
                    _ => {
                        return Err(ParseError::InvalidDefinitionParams(
                            "Cannot extend a definition that does not have transformations"
                                .to_string(),
                        ))
                    }
                };

                // add them to the list of transformations
                transformations.append(&mut original_transformations);
            } else {
                // this defines a new transformation
                let transformation = parse_transformation(item)?;
                transformations.push(transformation);
            }
        }

        Ok(Definition {
            name,
            value: DefinitionValue::Transformations(transformations),
        })
    } else {
        Err(ParseError::InvalidDefinitionParams(
            "Definition 'value' must be a mapping or a sequence".to_string(),
        ))
    }

    // parse the transformations
    // first, we have to get the 'value' on this definition
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

    let from = Point::new(from_values[0], from_values[1], from_values[2]);
    let to = Point::new(to_values[0], to_values[1], to_values[2]);
    let up = Vector::new(up_values[0], up_values[1], up_values[2]);

    let transformation = compute_view_transformation(from, to, up);
    let camera = Camera::new(width, height, field_of_view, transformation);

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

    let position = Point::new(position_values[0], position_values[1], position_values[2]);
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
    let reflectivity = get_or_default_f64!(description, "reflective", 0.);
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

fn construct_object(
    description: &Mapping,
    definitions: &[Definition],
    parent: Option<&group::Group>,
) -> Result<Arc<dyn Shape>> {
    let object_type = description
        .get("add")
        .ok_or(ParseError::InvalidValue(
            "An object definition must have an 'add' key".to_string(),
        ))?
        .as_str()
        .ok_or(ParseError::InvalidValue(
            "Object's 'add' must be a string".to_string(),
        ))?;

    let material = if let Some(material_description) = description.get("material") {
        // if the material is a string, it refers to a definition
        if let Some(definition_name) = material_description.as_str() {
            let definition = definitions
                .iter()
                .find(|x| x.name == definition_name)
                .ok_or(ParseError::InvalidValue(format!(
                    "No such defintion: {definition_name}"
                )))?;

            if let DefinitionValue::Material(material) = &definition.value {
                material.to_owned()
            } else {
                return Err(ParseError::InvalidValue(format!(
                    "Definition '{definition_name}' does not refer to a material"
                )));
            }
        } else {
            // otherwise, it's a material description
            parse_material(material_description)
        }
    } else {
        Material::default()
    };

    let mut transformation = Matrix::identity();

    if let Some(transformations) = description.get("transform").and_then(|x| x.as_sequence()) {
        for value in transformations {
            // values specified here can refer to definitions, and extend them

            if let Some(definition_name) = value.as_str() {
                let defined_transformations = &definitions
                    .iter()
                    .find(|x| x.name == definition_name)
                    .ok_or(ParseError::InvalidValue(format!(
                        "No such defintion: {definition_name}"
                    )))?
                    .value;

                if let DefinitionValue::Transformations(transformations) = defined_transformations {
                    for t in transformations {
                        match t {
                            Transformation::Scale(x, y, z) => {
                                transformation = transformation.scale(*x, *y, *z);
                            }
                            Transformation::Translate(x, y, z) => {
                                transformation = transformation.translate(*x, *y, *z);
                            }
                            Transformation::RotateX(radians) => {
                                transformation = transformation.rotate_x(*radians);
                            }
                            Transformation::RotateY(radians) => {
                                transformation = transformation.rotate_y(*radians);
                            }
                            Transformation::RotateZ(radians) => {
                                transformation = transformation.rotate_z(*radians);
                            }
                        }
                    }
                } else {
                    return Err(ParseError::InvalidValue(format!(
                        "Definition '{definition_name}' does not refer to a set of transformations"
                    )));
                }
            } else {
                // the transformations are specified directly
                let t = parse_transformation(value)?;
                match t {
                    Transformation::Scale(x, y, z) => {
                        transformation = transformation.scale(x, y, z);
                    }
                    Transformation::Translate(x, y, z) => {
                        transformation = transformation.translate(x, y, z);
                    }
                    Transformation::RotateX(radians) => {
                        transformation = transformation.rotate_x(radians);
                    }
                    Transformation::RotateY(radians) => {
                        transformation = transformation.rotate_y(radians);
                    }
                    Transformation::RotateZ(radians) => {
                        transformation = transformation.rotate_z(radians);
                    }
                }
            }
        }
    }

    let object: Arc<dyn Shape> = match object_type {
        "sphere" => Arc::new(sphere::Sphere::new(transformation, material, None)),
        "plane" => Arc::new(plane::Plane::new(transformation, material, None)),
        "cube" => Arc::new(cube::Cube::new(transformation, material, None)),
        "cylinder" => {
            // cylinders can be truncated and have a minimum and a maximum
            let minimum = get_or_default_f64!(description, "minimum", -f64::INFINITY);
            let maximum = get_or_default_f64!(description, "maximum", f64::INFINITY);
            let closed = description
                .get("closed")
                .and_then(|x| x.as_bool())
                .unwrap_or(false);

            Arc::new(cylinder::Cylinder::new(
                transformation,
                material,
                minimum,
                maximum,
                closed,
                None,
            ))
        }
        "cone" => {
            // ditto for cones
            let minimum = get_or_default_f64!(description, "minimum", -f64::INFINITY);
            let maximum = get_or_default_f64!(description, "maximum", f64::INFINITY);
            let closed = description
                .get("closed")
                .and_then(|x| x.as_bool())
                .unwrap_or(false);

            Arc::new(cone::Cone::new(
                transformation,
                material,
                minimum,
                maximum,
                closed,
                None,
            ))
        }
        "group" => {
            let group: Arc<dyn Shape> =
                Arc::new(group::Group::default().with_transformation(transformation));

            // if we have a parent given to this function
            // that means we are in a nested group.
            // We should be careful to add the group to the parent
            // before we add the children to the group,
            // because the group needs to be in the shapes store by the time we add children.
            // We must also have a unique reference to the group when we add it as a child
            // so we have to add this group to its parent at the level where it's defined.
            if let Some(parent) = parent {
                parent.add_child(&group);
            } else {
                // if we don't have a parent, then this is the root group
                // we need to just add it to the shapes store directly --
                // it won't be anyone's child
                let group_clone: Arc<dyn Shape> = Arc::clone(&group) as Arc<dyn Shape>;
                register_shape(group_clone);
            }

            let children = description
                .get("children")
                .ok_or(ParseError::InvalidGroupParams(
                    "Missing 'children' in group definition".to_string(),
                ))?
                .as_sequence()
                .ok_or(ParseError::InvalidGroupParams(
                    "Group 'children' must be a sequence".to_string(),
                ))?
                .iter()
                .map(|x| {
                    x.as_mapping().ok_or(ParseError::InvalidGroupParams(
                        "Each value in a group's 'children' must be a mapping".to_string(),
                    ))
                })
                .map(|x| construct_object(x?, definitions, group.as_group()))
                .collect::<Result<Vec<Arc<dyn Shape>>>>()?;

            // if the child is a group, then it's already been added
            // in the next level of recursion
            for child in children.iter() {
                if child.parent().is_none() {
                    let group = group.as_group().unwrap();
                    group.add_child(child);
                }
            }

            group
        }
        _ => unimplemented!(),
    };

    Ok(object)
}

pub fn parse_scene(input: &str) -> Result<Scene> {
    let sequence = serde_yaml::from_str::<Sequence>(input)?;

    let mut camera: Option<Camera> = None;
    let mut lights = vec![];
    let mut definitions = vec![];

    let mut objects: Vec<Arc<dyn Shape>> = vec![];

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
                    "sphere" | "plane" | "cube" | "cylinder" | "cone" | "group" => {
                        let object = construct_object(mapping, &definitions, None)?;
                        objects.push(object);
                    }
                    _ => return Err(ParseError::UnknownItem(item_type.to_owned())),
                }
            } else if let Some(Value::String(_)) = mapping.get("define") {
                let definition = parse_definition(mapping, &definitions)?;
                definitions.push(definition);
            }
        } else {
            return Err(ParseError::InvalidYamlObject(format!(
                "Expected mapping, found {item:#?}"
            )));
        }
    }

    let camera = camera.ok_or(ParseError::MissingCamera)?;

    let scene = Scene {
        camera,
        lights,
        objects,
    };
    Ok(scene)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_scene() {
        let input = r#"
        - add: camera
          width: 1000
          height: 500
          field-of-view: 1.0471975512
          from: [0, 5, -7]
          to: [0, 2, 0]
          up: [0, 1, 0]
        - add: light
          at: [-10, 10, -10]
          intensity: [1, 1, 1]
        - add: light
          at: [10, 10, -10]
          intensity: [1, 1, 1]
        - add: sphere
          material:
            color: [1, 0.2, 1]
            diffuse: 0.7
            specular: 0.3
          transform:
            - [ scale, 1, 0.5, 1 ]
            - [translate, 0, 1, 0]
        - add: plane
          material:
            color: [1, 0.9, 0.9]
            specular: 0
          transform:
            - [rotate-x, 1.5708]
            - [ translate, 0, 0, 5]
        "#;

        let scene = parse_scene(input).unwrap();

        assert_eq!(scene.lights.len(), 2);
        assert_eq!(scene.lights[0].position, Point::new(-10., 10., -10.));
        assert_eq!(scene.lights[0].intensity, Color::new(1., 1., 1.));
        assert_eq!(scene.lights[1].position, Point::new(10., 10., -10.));
        assert_eq!(scene.lights[1].intensity, Color::new(1., 1., 1.));

        assert_eq!(scene.objects.len(), 2);
        assert_eq!(scene.objects[0].material().color, Color::new(1., 0.2, 1.));
        assert_eq!(scene.objects[0].material().diffuse, 0.7);
        assert_eq!(scene.objects[0].material().specular, 0.3);
        assert_eq!(
            scene.objects[0].transformation(),
            Matrix::identity().scale(1., 0.5, 1.).translate(0., 1., 0.)
        );
    }

    #[test]
    fn parse_scene_with_definitions() {
        let input = r#"
        - add: camera
          width: 1000
          height: 500
          field-of-view: 1.0471975512
          from: [0, 5, -7]
          to: [0, 2, 0]
          up: [0, 1, 0]
        - add: light
          at: [-10, 10, -10]
          intensity: [1, 1, 1]
        - define: white-material
          value:
            color: [1, 1, 1]
            diffuse: 0.7
            ambient: 0.1
            shininess: 250
            specular: 0.0
            reflective: 0.1
        - define: blue-material
          extend: white-material
          value:
            color: [0.537, 0.831, 0.914]
        - define: standard-transform
          value:
            - [ translate, 1, -1, 1 ]
            - [ scale, 0.5, 0.5, 0.5 ]
        - define: large-object
          value:
            - standard-transform
            - [ scale, 3.5, 3.5, 3.5 ]
        - add: cube
          material: blue-material
          transform:
            - large-object
            - [ translate, 0.5, 1.5, -0.5 ]
        "#;

        let scene = parse_scene(input).unwrap();
        assert_eq!(scene.objects.len(), 1);
        assert_eq!(
            scene.objects[0].material().color,
            Color::new(0.537, 0.831, 0.914)
        );

        assert_eq!(scene.objects[0].material().diffuse, 0.7);
        assert_eq!(scene.objects[0].material().ambient, 0.1);
        assert_eq!(scene.objects[0].material().specular, 0.0);
        assert_eq!(scene.objects[0].material().shininess, 250.);
        assert_eq!(scene.objects[0].material().reflectivity, 0.1);
        assert_eq!(
            scene.objects[0].transformation(),
            Matrix::identity()
                .translate(1., -1., 1.)
                .scale(0.5, 0.5, 0.5)
                .scale(3.5, 3.5, 3.5)
                .translate(0.5, 1.5, -0.5)
        );
    }

    #[test]
    fn parse_scene_with_group() {
        let input = r#"
        - add: camera
          width: 1000
          height: 500
          field-of-view: 1.0471975512
          from: [0, 5, -7]
          to: [0, 2, 0]
          up: [0, 1, 0]
        - add: light
          at: [-10, 10, -10]
          intensity: [1, 1, 1]
        - add: group
          transform:
            - [ translate, 0, 1, 0 ]
          children:
            - add: sphere
              material:
                color: [1, 0.2, 1]
                diffuse: 0.7
                specular: 0.3
              transform:
                - [ scale, 1, 0.5, 1 ]
            - add: plane
              material:
                color: [1, 0.9, 0.9]
                specular: 0
              transform:
                - [rotate-x, 1.5708]
                - [ translate, 0, 0, 5]
        "#;

        let scene = parse_scene(input).unwrap();
        assert_eq!(scene.objects.len(), 1);

        let group = scene.objects[0].as_group().unwrap();
        assert_eq!(
            group.transformation(),
            Matrix::identity().translate(0., 1., 0.)
        );

        let sphere = group.get_child(0).unwrap();
        assert_eq!(
            sphere.material(),
            Material {
                color: Color::new(1., 0.2, 1.),
                diffuse: 0.7,
                specular: 0.3,
                ..Default::default()
            }
        );
        assert_eq!(sphere.transformation(), Matrix::scaling(1., 0.5, 1.));

        let plane = group.get_child(1).unwrap();
        assert_eq!(
            plane.material(),
            Material {
                color: Color::new(1., 0.9, 0.9),
                specular: 0.,
                ..Default::default()
            }
        );
        assert_eq!(
            plane.transformation(),
            Matrix::identity().rotate_x(1.5708).translate(0., 0., 5.)
        );
    }

    #[test]
    fn parse_scene_with_nested_groups() {
        let input = r#"
        - add: camera
          width: 1000
          height: 500
          field-of-view: 1.0471975512
          from: [0, 5, -7]
          to: [0, 2, 0]
          up: [0, 1, 0]
        - add: light
          at: [-10, 10, -10]
          intensity: [1, 1, 1]
        - add: group
          transform:
            - [ translate, 0, 1, 0 ]
          children:
            - add: group
              transform:
                - [ scale, 1, 0.5, 1 ]
              children:
                - add: sphere
                  material:
                    color: [1, 0.2, 1]
                    diffuse: 0.7
                    specular: 0.3
                - add: cube
                  material: {}
            - add: plane
              material:
                color: [1, 0.9, 0.9]
                specular: 0
              transform:
                - [rotate-x, 1.5708]
                - [ translate, 0, 0, 5]
        "#;

        let scene = parse_scene(input).unwrap();
        assert_eq!(scene.objects.len(), 1);
        assert_eq!(
            scene.objects[0].transformation(),
            Matrix::identity().translate(0., 1., 0.)
        );

        let group = scene.objects[0].as_group().unwrap();

        let child_group = group.get_child(0).unwrap();
        let child_group = child_group.as_group().unwrap();
        let sphere = child_group.get_child(0).unwrap();
        assert_eq!(
            sphere.material(),
            Material {
                color: Color::new(1., 0.2, 1.),
                diffuse: 0.7,
                specular: 0.3,
                ..Default::default()
            }
        );

        let cube = child_group.get_child(1).unwrap();
        assert_eq!(cube.material(), Material::default());

        let plane = group.get_child(1).unwrap();
        assert_eq!(
            plane.material(),
            Material {
                color: Color::new(1., 0.9, 0.9),
                specular: 0.,
                ..Default::default()
            }
        );
        assert_eq!(
            plane.transformation(),
            Matrix::identity().rotate_x(1.5708).translate(0., 0., 5.)
        );
    }
}
