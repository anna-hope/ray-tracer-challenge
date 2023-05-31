use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::sync::Arc;

use thiserror::Error;

use ray_tracer::prelude::*;
use ray_tracer::shape::{group::Group, smooth_triangle::SmoothTriangle, triangle::Triangle};

#[derive(Debug, Error)]
pub enum ParserError {
    #[error("failed to open the file: {0}")]
    CouldNotOpenFile(#[from] std::io::Error),

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("missing vertex at index {0}")]
    MissingVertex(usize),
}

pub type Result<T> = std::result::Result<T, ParserError>;

#[derive(Debug, Clone)]
pub struct ParsedObj {
    pub ignored_lines: u32,
    pub default_group: Arc<Group>,
    pub vertices: Vec<Point>,
    pub normals: Vec<Vector>,
    groups: HashMap<String, Arc<dyn Shape>>,
}

impl ParsedObj {
    pub fn get_group(&self, group_name: &str) -> Option<&Group> {
        self.groups.get(group_name)?.as_group()
    }
}

fn parse_obj_string(value: impl Into<String>, material: Option<Material>) -> Result<ParsedObj> {
    let string: String = value.into();

    let mut ignored_lines = 0;
    let mut vertices = vec![];
    let mut normals = vec![];

    let default_group =
        Arc::new(Group::default().with_material(material.clone().unwrap_or_default()));
    let default_group_clone = Arc::clone(&default_group) as Arc<dyn Shape>;
    register_shape(default_group_clone);

    let mut groups = HashMap::new();
    let mut current_group: Option<Arc<dyn Shape>> = None;

    for (n, line) in string.lines().enumerate() {
        let tokens = line.split_whitespace().collect::<Vec<_>>();
        if tokens.is_empty() {
            continue;
        }

        let line_no = n + 1;

        match tokens[0] {
            "v" | "vn" => {
                if let (Some(x), Some(y), Some(z)) = (tokens.get(1), tokens.get(2), tokens.get(3)) {
                    let x = x.parse::<f64>().map_err(|_| {
                        ParserError::InvalidInput(format!(
                            "Invalid vertex record on line {line_no}: {x}"
                        ))
                    })?;
                    let y = y.parse::<f64>().map_err(|_| {
                        ParserError::InvalidInput(format!(
                            "Invalid vertex record on line {line_no}: {y}"
                        ))
                    })?;
                    let z = z.parse::<f64>().map_err(|_| {
                        ParserError::InvalidInput(format!(
                            "Invalid vertex record on line {line_no}: {z}"
                        ))
                    })?;

                    if tokens[0] == "v" {
                        vertices.push(Point::new(x, y, z));
                    } else {
                        // vn
                        normals.push(Vector::new(x, y, z));
                    }
                } else {
                    println!("Malformed vertex record on line {line_no}");
                }
            }
            "f" => {
                let mut face_vertices = vec![];
                let mut normal_indices = vec![];

                for token in &tokens[1..] {
                    let subtokens = token.split('/').collect::<Vec<_>>();

                    let vertex_index = subtokens.first().ok_or(ParserError::InvalidInput(
                        format!("Missing vertex index on line {line_no}"),
                    ))?;
                    let vertex_index = vertex_index.parse::<usize>().map_err(|_| {
                        ParserError::InvalidInput(format!(
                            "Invalid vertex index for face record on line {n}: {vertex_index}"
                        ))
                    })?;

                    // subtract 1 because indices in OBJ files are 1-based
                    let vertex_index = vertex_index - 1;
                    let vertex = vertices
                        .get(vertex_index)
                        .ok_or(ParserError::MissingVertex(vertex_index + 1))?;
                    face_vertices.push(*vertex);

                    if let Some(normal_index) = subtokens.get(2) {
                        let normal_index = normal_index.parse::<usize>().map_err(|_| {
                            ParserError::InvalidInput(format!(
                                "Invalid normal index for face record on line {n}: {normal_index}"
                            ))
                        })?;
                        normal_indices.push(normal_index - 1);
                    }
                }

                let mut triangles = vec![];
                for index in 1..face_vertices.len() - 1 {
                    let point1 = face_vertices[0];
                    let point2 = face_vertices[index];
                    let point3 = face_vertices[index + 1];

                    let triangle = if normal_indices.is_empty() {
                        let triangle = Triangle::new(point1, point2, point3, Material::default());
                        Arc::new(triangle) as Arc<dyn Shape>
                    } else {
                        let normal1 = normals[normal_indices[0]];
                        let normal2 = normals[normal_indices[1]];
                        let normal3 = normals[normal_indices[2]];
                        let triangle = SmoothTriangle::new(
                            point1,
                            point2,
                            point3,
                            normal1,
                            normal2,
                            normal3,
                            Material::default(),
                        );
                        Arc::new(triangle) as Arc<dyn Shape>
                    };

                    triangles.push(triangle);
                }

                if let Some(ref group) = current_group {
                    let group = group.as_group().unwrap();

                    for triangle in triangles.iter() {
                        group.add_child(triangle);
                    }
                } else {
                    for triangle in triangles.iter() {
                        default_group.add_child(triangle);
                    }
                }
            }
            "g" => {
                let group_name = tokens
                    .get(1)
                    .ok_or(ParserError::InvalidInput(format!(
                        "Malformed named group record on line {n}: missing group name"
                    )))?
                    .to_string();

                let group: Arc<dyn Shape> =
                    Arc::new(Group::default().with_material(material.clone().unwrap_or_default()));
                default_group.add_child(&group);
                groups.insert(group_name, Arc::clone(&group));
                current_group = Some(group);
            }
            _ => {
                ignored_lines += 1;
            }
        }
    }

    Ok(ParsedObj {
        ignored_lines,
        vertices,
        normals,
        default_group,
        groups,
    })
}

pub fn parse_obj_file(filename: &str, material: Option<Material>) -> Result<ParsedObj> {
    let mut file = File::open(filename)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;

    parse_obj_string(buffer, material)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ignore_unrecognized_statements() {
        let gibberish = r#"There was a young lady named Bright
        Who traveled much faster than light.
        She set out one day
        In a relative way,
        And came back the previous night."#;

        let parsed = parse_obj_string(gibberish, None).unwrap();
        assert_eq!(parsed.ignored_lines, 5);
    }

    #[test]
    fn process_vertex_data() {
        let obj_string = r#"
#
# some random gibberish here
#

# note the two spaces in the line below

v  -1 1 0
v -1.0000 0.5000 0.0000
v 1 0 0
v 1 1 0
        "#;

        let parsed = parse_obj_string(obj_string, None).unwrap();
        assert_eq!(parsed.vertices[0], Point::new(-1., 1., 0.));
        assert_eq!(parsed.vertices[1], Point::new(-1., 0.5, 0.));
        assert_eq!(parsed.vertices[2], Point::new(1., 0., 0.));
        assert_eq!(parsed.vertices[3], Point::new(1., 1., 0.));
    }

    #[test]
    fn parse_triangle_faces() {
        let obj_string = r#"
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0

f 1 2 3
f 1 3 4
            "#;

        let parsed = parse_obj_string(obj_string, None).unwrap();
        let group = parsed.default_group;
        let child1 = group.get_child(0).unwrap();
        let triangle1 = child1.as_any().downcast_ref::<Triangle>().unwrap();

        let child2 = group.get_child(1).unwrap();

        let triangle2 = child2.as_any().downcast_ref::<Triangle>().unwrap();

        assert_eq!(triangle1.point1, parsed.vertices[0]);
        assert_eq!(triangle1.point2, parsed.vertices[1]);
        assert_eq!(triangle1.point3, parsed.vertices[2]);

        assert_eq!(triangle2.point1, parsed.vertices[0]);
        assert_eq!(triangle2.point2, parsed.vertices[2]);
        assert_eq!(triangle2.point3, parsed.vertices[3]);
    }

    #[test]
    fn triangulate_polygons() {
        let obj_string = r#"
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0
v 0 2 0

f 1 2 3 4 5
        "#;

        let parsed = parse_obj_string(obj_string, None).unwrap();
        let group = parsed.default_group;

        let child1 = group.get_child(0).unwrap();
        let triangle1 = child1.as_any().downcast_ref::<Triangle>().unwrap();

        let child2 = group.get_child(1).unwrap();
        let triangle2 = child2.as_any().downcast_ref::<Triangle>().unwrap();

        let child3 = group.get_child(2).unwrap();
        let triangle3 = child3.as_any().downcast_ref::<Triangle>().unwrap();

        assert_eq!(triangle1.point1, parsed.vertices[0]);
        assert_eq!(triangle1.point2, parsed.vertices[1]);
        assert_eq!(triangle1.point3, parsed.vertices[2]);

        assert_eq!(triangle2.point1, parsed.vertices[0]);
        assert_eq!(triangle2.point2, parsed.vertices[2]);
        assert_eq!(triangle2.point3, parsed.vertices[3]);

        assert_eq!(triangle3.point1, parsed.vertices[0]);
        assert_eq!(triangle3.point2, parsed.vertices[3]);
        assert_eq!(triangle3.point3, parsed.vertices[4]);
    }

    #[test]
    fn triangles_in_groups() {
        let filename = "triangles.obj";
        let parsed = parse_obj_file(filename, None).unwrap();
        let group1 = parsed.get_group("FirstGroup").unwrap();
        let group2 = parsed.get_group("SecondGroup").unwrap();

        let group1_child1 = group1.get_child(0).unwrap();
        let triangle1 = group1_child1.as_any().downcast_ref::<Triangle>().unwrap();

        let group2_child1 = group2.get_child(0).unwrap();
        let triangle2 = group2_child1.as_any().downcast_ref::<Triangle>().unwrap();

        assert_eq!(triangle1.point1, parsed.vertices[0]);
        assert_eq!(triangle1.point2, parsed.vertices[1]);
        assert_eq!(triangle1.point3, parsed.vertices[2]);

        assert_eq!(triangle2.point1, parsed.vertices[0]);
        assert_eq!(triangle2.point2, parsed.vertices[2]);
        assert_eq!(triangle2.point3, parsed.vertices[3]);
    }

    #[test]
    fn convert_obj_file_to_group() {
        let filename = "triangles.obj";
        let parsed = parse_obj_file(filename, None).unwrap();
        let group = Arc::clone(&parsed.default_group);

        let child1 = group.get_child(0).unwrap();
        let child_group1 = child1.as_group().unwrap();

        let child2 = group.get_child(1).unwrap();
        let child_group2 = child2.as_group().unwrap();

        let first_group = parsed.get_group("FirstGroup").unwrap();
        let second_group = parsed.get_group("SecondGroup").unwrap();

        assert_eq!(child_group1, first_group);
        assert_eq!(child_group2, second_group);
    }

    #[test]
    fn vertex_normal_records() {
        let obj_string = r#"
vn 0 0 1
vn 0.707 0 -0.707
vn 1 2 3
"#;
        let parsed = parse_obj_string(obj_string, None).unwrap();
        assert_eq!(parsed.normals[0], Vector::new(0., 0., 1.));
        assert_eq!(parsed.normals[1], Vector::new(0.707, 0., -0.707));
        assert_eq!(parsed.normals[2], Vector::new(1., 2., 3.));
    }

    #[test]
    fn faces_with_normals() {
        let obj_string = r#"
v 0 1 0
v -1 0 0
v 1 0 0

vn -1 0 0
vn 1 0 0
vn 0 1 0

f 1//3 2//1 3//2
f 1/0/3 2/102/1 3/14/2
"#;

        let parsed = parse_obj_string(obj_string, None).unwrap();
        let group = parsed.default_group;

        let group1_child1 = group.get_child(0).unwrap();
        let triangle1 = group1_child1
            .as_any()
            .downcast_ref::<SmoothTriangle>()
            .unwrap();

        let group2_child1 = group.get_child(1).unwrap();
        let triangle2 = group2_child1
            .as_any()
            .downcast_ref::<SmoothTriangle>()
            .unwrap();

        assert_eq!(triangle1.point1, parsed.vertices[0]);
        assert_eq!(triangle1.point2, parsed.vertices[1]);
        assert_eq!(triangle1.point3, parsed.vertices[2]);

        assert_eq!(triangle1.normal1, parsed.normals[2]);
        assert_eq!(triangle1.normal2, parsed.normals[0]);
        assert_eq!(triangle1.normal3, parsed.normals[1]);

        assert_eq!(triangle2.point1, parsed.vertices[0]);
        assert_eq!(triangle2.point2, parsed.vertices[1]);
        assert_eq!(triangle2.point3, parsed.vertices[2]);

        assert_eq!(triangle2.normal1, parsed.normals[2]);
        assert_eq!(triangle2.normal2, parsed.normals[0]);
        assert_eq!(triangle2.normal3, parsed.normals[1]);
    }
}
