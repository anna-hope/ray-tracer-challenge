use std::fmt::Debug;

use crate::{shape::Shape, Matrix, Result, Tuple, EPSILON};

/// Ray.origin is a point, Ray.direction is a vector.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Tuple,
    pub direction: Tuple,
}

impl Ray {
    /// Origin should be a point, direction should be a vector.
    pub fn new(origin: Tuple, direction: Tuple) -> Self {
        Self { origin, direction }
    }

    pub fn position(&self, t: f64) -> Tuple {
        self.origin + self.direction * t
    }

    pub fn transform(&self, matrix: Matrix) -> Self {
        let origin = matrix.clone() * self.origin;
        let direction = matrix * self.direction;
        Self { origin, direction }
    }
}

pub trait Intersect {
    /// Computes the intersections of this ray with the given object.
    fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>>;
}

#[derive(Debug)]
pub struct Computations {
    pub t: f64,
    pub object: Box<dyn Shape>,
    pub point: Tuple,
    pub eye_vector: Tuple,
    pub normal_vector: Tuple,
    pub inside: bool,
    pub over_point: Tuple,
    pub reflect_vector: Tuple,
    pub n1: Option<f64>,
    pub n2: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct Intersection {
    pub t: f64,
    pub object: Box<dyn Shape>,
}

impl Intersection {
    pub fn new(t: f64, object: Box<dyn Shape>) -> Self {
        Self { t, object }
    }

    /// Precomputes the point (in world space) where the intersection occurred
    /// the eye vector (pointing back toward the eye/camera)
    /// and the normal vector.
    pub fn prepare_computations(
        &self,
        ray: &Ray,
        intersections: &[Intersection],
    ) -> Result<Computations> {
        let point = ray.position(self.t);
        let eye_vector = -ray.direction;
        let normal_vector = self.object.normal_at(point)?;

        // if the dot product of normal_vector and eye_vector is negative
        // then they're pointing in (roughly) opposite directions
        let normal_dot_eye = normal_vector.dot(&eye_vector)?;
        let (inside, normal_vector) = if normal_dot_eye < 0. {
            (true, -normal_vector)
        } else {
            (false, normal_vector)
        };

        let reflect_vector = ray.direction.reflect(normal_vector)?;
        let over_point = point + normal_vector * EPSILON;

        // I'm sorry for this
        let mut containers: Vec<&Box<dyn Shape>> = vec![];
        let mut n1: Option<f64> = None;
        let mut n2: Option<f64> = None;

        if let Some(hit) = hit(&intersections.clone()) {
            for intersection in intersections {
                if intersection == hit {
                    if containers.is_empty() {
                        n1 = Some(1.);
                    } else {
                        if let Some(object) = containers.last() {
                            n1 = Some(object.material().refractive_index);
                        }
                    }
                }

                let object = &intersection.object;
                if let Some(position) = containers.iter().position(|&x| x == object) {
                    containers.remove(position);
                } else {
                    containers.push(object);
                }

                if intersection == hit {
                    if containers.is_empty() {
                        n2 = Some(1.);
                    } else {
                        if let Some(object) = containers.last() {
                            n2 = Some(object.material().refractive_index);
                        }
                    }

                    break;
                }
            }
        }

        Ok(Computations {
            t: self.t,
            object: self.object.clone(),
            point,
            eye_vector,
            normal_vector,
            inside,
            over_point,
            reflect_vector,
            n1,
            n2,
        })
    }
}

impl PartialEq for Intersection {
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t
            && self.object.id() == other.object.id()
            && self.object.shape_type() == other.object.shape_type()
    }
}

// impl<'a> Debug for Intersection<'a> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         f.debug_struct("Intersection")
//             .field("t", &self.t)
//             .field("object type", &self.object.shape_type())
//             .field("object id", &self.object.id())
//             .finish()
//     }
// }

/// Finds the intersection that hits the object.
/// Always picks the lowest non-negative intersection (intersection with the smallest t > 0).
pub fn hit(intersections: &[Intersection]) -> Option<&Intersection> {
    let mut lowest_nonnegative_intersection: Option<&Intersection> = None;
    for intersection in intersections {
        if let Some(i) = lowest_nonnegative_intersection {
            if intersection.t > 0. && intersection.t < i.t {
                lowest_nonnegative_intersection = Some(intersection);
            }
        } else if intersection.t > 0. {
            lowest_nonnegative_intersection = Some(intersection);
        }
    }
    lowest_nonnegative_intersection
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        shape::{plane::Plane, sphere::Sphere},
        EPSILON,
    };

    #[test]
    fn create_and_query_ray() {
        let origin = Tuple::point(1., 2., 3.);
        let direction = Tuple::vector(4., 5., 6.);
        let ray = Ray::new(origin, direction);
        assert_eq!(ray.origin, origin);
        assert_eq!(ray.direction, direction);
    }

    #[test]
    fn compute_point_from_distance() {
        let ray = Ray::new(Tuple::point(2., 3., 4.), Tuple::vector(1., 0., 0.));
        assert_eq!(ray.position(0.), Tuple::point(2., 3., 4.));
        assert_eq!(ray.position(1.), Tuple::point(3., 3., 4.));
        assert_eq!(ray.position(-1.), Tuple::point(1., 3., 4.));
        assert_eq!(ray.position(2.5), Tuple::point(4.5, 3., 4.));
    }

    #[test]
    fn intersection_encapsulates_t_and_object() {
        let sphere = Sphere::default();
        let intersection = Intersection::new(3.5, Box::new(sphere.clone()));
        assert_eq!(intersection.t, 3.5);
        // assert_eq!(&sphere, intersection.object);
        assert_eq!(intersection.object.id(), sphere.id());
        assert_eq!(intersection.object.shape_type(), sphere.shape_type());
    }

    #[test]
    fn aggregating_intersections() {
        let sphere = Sphere::default();
        let i1 = Intersection::new(1., Box::new(sphere.clone()));
        let i2 = Intersection::new(2., Box::new(sphere));
        let xs = &[i1, i2];
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_all_intersections_have_positive_t() {
        let sphere = Sphere::default();
        let i1 = Intersection::new(1., Box::new(sphere.clone()));
        let i2 = Intersection::new(2., Box::new(sphere));
        let xs = &[i1.clone(), i2.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i1);
    }

    #[test]
    fn hit_some_intersections_have_negative_t() {
        let sphere = Sphere::default();
        let i1 = Intersection::new(-1., Box::new(sphere.clone()));
        let i2 = Intersection::new(1., Box::new(sphere.clone()));
        let xs = &[i1.clone(), i2.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i2);
    }

    #[test]
    fn hit_all_intersections_have_negative_t() {
        let sphere = Box::new(Sphere::default());
        let i1 = Intersection::new(-2., sphere.clone());
        let i2 = Intersection::new(-1., sphere);
        let xs = &[i1, i2];
        let intersection = hit(xs);
        assert!(intersection.is_none());
    }

    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let sphere = Box::new(Sphere::default());
        let i1 = Intersection::new(5., sphere.clone());
        let i2 = Intersection::new(7., sphere.clone());
        let i3 = Intersection::new(-3., sphere.clone());
        let i4 = Intersection::new(2., sphere);
        let xs = &[i1.clone(), i2.clone(), i3.clone(), i4.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i4);
    }

    #[test]
    fn translate_ray() {
        let ray = Ray::new(Tuple::point(1., 2., 3.), Tuple::vector(0., 1., 0.));
        let matrix = Matrix::translation(3., 4., 5.);
        let ray2 = ray.transform(matrix);
        assert_eq!(ray2.origin, Tuple::point(4., 6., 8.));
        assert_eq!(ray2.direction, Tuple::vector(0., 1., 0.));
    }

    #[test]
    fn scale_ray() {
        let ray = Ray::new(Tuple::point(1., 2., 3.), Tuple::vector(0., 1., 0.));
        let matrix = Matrix::scaling(2., 3., 4.);
        let ray2 = ray.transform(matrix);
        assert_eq!(ray2.origin, Tuple::point(2., 6., 12.));
        assert_eq!(ray2.direction, Tuple::vector(0., 3., 0.));
    }

    #[test]
    fn precompute_state_of_intersection() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = Box::new(Sphere::default());
        let intersection = Intersection::new(4., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert_eq!(comps.t, intersection.t);
        assert_eq!(comps.object.id(), intersection.object.id());
        assert_eq!(comps.object.shape_type(), intersection.object.shape_type());
        assert_eq!(comps.point, Tuple::point(0., 0., -1.));
        assert_eq!(comps.eye_vector, Tuple::vector(0., 0., -1.));
        assert_eq!(comps.normal_vector, Tuple::vector(0., 0., -1.));
    }

    #[test]
    fn hit_when_intersection_outside() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = Box::new(Sphere::default());
        let intersection = Intersection::new(4., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert!(!comps.inside);
    }

    #[test]
    fn hit_when_intersection_inside() {
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let shape = Box::new(Sphere::default());
        let intersection = Intersection::new(1., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert_eq!(comps.point, Tuple::point(0., 0., 1.));
        assert_eq!(comps.eye_vector, Tuple::vector(0., 0., -1.));
        assert!(comps.inside);
        assert_eq!(comps.normal_vector, Tuple::vector(0., 0., -1.));
    }

    #[test]
    fn hit_should_offset_point() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape =
            Box::new(Sphere::default().with_transformation(Matrix::translation(0., 0., 1.)));
        let intersection = Intersection::new(5., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert!(comps.over_point.z < -EPSILON / 2.);
        assert!(comps.point.z > comps.over_point.z);
    }

    #[test]
    fn precompute_reflection_vector() {
        let shape = Box::new(Plane::new());
        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Tuple::point(0., 1., -1.), Tuple::vector(0., -val, val));
        let intersection = Intersection::new(2.0_f64.sqrt(), shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert_eq!(comps.reflect_vector, Tuple::vector(0., val, val));
    }

    #[test]
    fn find_n1_and_n2_at_various_intersections() {
        let mut a = Box::new(Sphere::glass().with_transformation(Matrix::scaling(2., 2., 2.)));
        let mut a_material = a.material();
        a_material.refractive_index = 1.5;
        a.set_material(a_material);

        let mut b =
            Box::new(Sphere::glass().with_transformation(Matrix::translation(0., 0., -0.25)));
        let mut b_material = b.material();
        b_material.refractive_index = 2.;
        b.set_material(b_material);

        let mut c =
            Box::new(Sphere::glass().with_transformation(Matrix::translation(0., 0., 0.25)));
        let mut c_material = c.material();
        c_material.refractive_index = 2.5;
        c.set_material(c_material);

        let ray = Ray::new(Tuple::point(0., 0., -4.), Tuple::vector(0., 0., 1.));
        let xs = vec![
            Intersection::new(2., a.clone()),
            Intersection::new(2.75, b.clone()),
            Intersection::new(3.25, c.clone()),
            Intersection::new(4.75, b),
            Intersection::new(5.25, c),
            Intersection::new(6., a),
        ];

        let n1s = [1., 1.5, 2., 2.5, 2.5, 1.5];
        let n2s = [1.5, 2., 2.5, 2.5, 1.5, 1.];

        for ((n1, n2), intersection) in n1s.into_iter().zip(n2s).zip(&xs) {
            let comps = intersection.prepare_computations(&ray, &xs).unwrap();
            assert_eq!(comps.n1, Some(n1));
            assert_eq!(comps.n2, Some(n2));
        }
    }
}
