use std::fs::File;
use std::io::Read;
use std::sync::Arc;

use thiserror::Error;

use ray_tracer::prelude::*;
use ray_tracer::shape::{group::Group, triangle::Triangle};

#[derive(Debug, Error)]
pub enum ParserError {
    #[error("failed to open the file")]
    CouldNotOpenFile(#[from] std::io::Error),

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("missing vertex at index {0}")]
    MissingVertex(usize),
}

pub type Result<T> = std::result::Result<T, ParserError>;

#[derive(Debug, Clone)]
pub struct ParsedObj {
    pub ignored_lines: usize,
    pub vertices: Vec<Point>,
    pub default_group: Arc<Group>,
}

fn apply_fan_triangulation(vertices: &[Point]) -> Vec<Triangle> {
    let mut triangles = vec![];

    for index in 1..vertices.len() - 1 {
        let triangle = Triangle::new(
            vertices[0],
            vertices[index],
            vertices[index + 1],
            Material::default(),
        );
        triangles.push(triangle);
    }

    triangles
}

fn parse_obj_string(value: impl Into<String>) -> Result<ParsedObj> {
    let string: String = value.into();

    let mut ignored_lines = 0;
    let mut vertices = vec![];

    let default_group = Arc::new(Group::default());
    let default_group_clone = Arc::clone(&default_group) as Arc<dyn Shape>;
    register_shape(default_group_clone);

    for (n, line) in string.lines().enumerate() {
        let items = line.trim().split(' ').collect::<Vec<_>>();
        if items.is_empty() {
            continue;
        }

        match items[0] {
            "v" => {
                if let (Some(x), Some(y), Some(z)) = (items.get(1), items.get(2), items.get(3)) {
                    let x = x.parse::<f64>().map_err(|_| {
                        ParserError::InvalidInput(format!("Invalid vertex record on line {n}: {x}"))
                    })?;
                    let y = y.parse::<f64>().map_err(|_| {
                        ParserError::InvalidInput(format!("Invalid vertex record on line {n}: {y}"))
                    })?;
                    let z = z.parse::<f64>().map_err(|_| {
                        ParserError::InvalidInput(format!("Invalid vertex record on line {n}: {z}"))
                    })?;

                    vertices.push(Point::new(x, y, z));
                } else {
                    println!("Malformed vertex record on line {n}");
                }
            }
            "f" => {
                let mut face_vertices = vec![];
                for vertex_index in &items[1..] {
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
                }

                let mut triangles = apply_fan_triangulation(&face_vertices)
                    .iter()
                    .map(|x| Arc::new(x.to_owned()) as Arc<dyn Shape>)
                    .collect::<Vec<_>>();

                for mut triangle in triangles.iter_mut() {
                    default_group.add_child(&mut triangle);
                }
            }
            "vn" => todo!(),
            _ => {
                ignored_lines += 1;
            }
        }
    }

    Ok(ParsedObj {
        ignored_lines,
        vertices,
        default_group,
    })
}

pub fn parse_obj_file(filename: &str) -> Result<ParsedObj> {
    let mut file = File::open(filename)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;

    parse_obj_string(buffer)
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

        let parsed = parse_obj_string(gibberish).unwrap();
        assert_eq!(parsed.ignored_lines, 5);
    }

    #[test]
    fn process_vertex_data() {
        let obj_string = r#"
v -1 1 0
v -1.0000 0.5000 0.0000
v 1 0 0
v 1 1 0
        "#;

        let parsed = parse_obj_string(obj_string).unwrap();
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

        let parsed = parse_obj_string(obj_string).unwrap();
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

        let parsed = parse_obj_string(obj_string).unwrap();
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
}
