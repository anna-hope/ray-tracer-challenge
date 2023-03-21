use std::fmt::Debug;

use crate::{Matrix, SceneObject, Tuple};

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
    fn intersect(&self, ray: &Ray) -> Vec<Intersection>;
}

#[derive(Clone)]
pub struct Intersection<'a> {
    pub t: f64,
    pub object: &'a dyn SceneObject,
}

impl<'a> Intersection<'a> {
    pub fn new(t: f64, object: &'a dyn SceneObject) -> Self {
        Self { t, object }
    }
}

impl<'a> PartialEq for Intersection<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t
            && self.object.id() == other.object.id()
            && self.object.object_type() == other.object.object_type()
    }
}

impl<'a> Debug for Intersection<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Intersection")
            .field("t", &self.t)
            .field("object type", &self.object.object_type())
            .field("object id", &self.object.id())
            .finish()
    }
}

/// Finds the intersection that hits the object.
/// Always picks the lowest non-negative intersection (intersection with the smallest t > 0).
pub fn hit<'a>(xs: &'a [Intersection]) -> Option<&'a Intersection<'a>> {
    let mut lowest_nonnegative_intersection: Option<&Intersection> = None;
    for intersection in xs {
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
    use crate::sphere::Sphere;

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
        let sphere = Sphere::new();
        let intersection = Intersection::new(3.5, &sphere);
        assert_eq!(intersection.t, 3.5);
        // assert_eq!(&sphere, intersection.object);
        assert_eq!(intersection.object.id(), sphere.id());
        assert_eq!(intersection.object.object_type(), sphere.object_type());
    }

    #[test]
    fn aggregating_intersections() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(1., &sphere);
        let i2 = Intersection::new(2., &sphere);
        let xs = &[i1.clone(), i2.clone()];
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_all_intersections_have_positive_t() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(1., &sphere);
        let i2 = Intersection::new(2., &sphere);
        let xs = &[i1.clone(), i2.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i1);
    }

    #[test]
    fn hit_some_intersections_have_negative_t() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(-1., &sphere);
        let i2 = Intersection::new(1., &sphere);
        let xs = &[i1.clone(), i2.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i2);
    }

    #[test]
    fn hit_all_intersections_have_negative_t() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(-2., &sphere);
        let i2 = Intersection::new(-1., &sphere);
        let xs = &[i1, i2];
        let intersection = hit(xs);
        assert!(intersection.is_none());
    }

    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(5., &sphere);
        let i2 = Intersection::new(7., &sphere);
        let i3 = Intersection::new(-3., &sphere);
        let i4 = Intersection::new(2., &sphere);
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
}
