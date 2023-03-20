use crate::Tuple;

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
}

pub trait Intersect<T> {
    fn intersect(&self, ray: &Ray) -> Vec<Intersection<T>>;
}

#[derive(Debug, PartialEq, PartialOrd)]
pub struct Intersection<'a, T> {
    pub t: f64,
    pub object: &'a T,
}

impl<'a, T> Intersection<'a, T> {
    pub fn new(t: f64, object: &'a T) -> Self {
        Self { t, object }
    }
}

pub fn hit<'a, T>(xs: &'a [&'a Intersection<T>]) -> Option<&'a Intersection<'a, T>> {
    let mut lowest_nonnegative_intersection: Option<&Intersection<T>> = None;
    for intersection in xs {
        if let Some(i) = lowest_nonnegative_intersection {
            if intersection.t > 0. && intersection.t < i.t {
                lowest_nonnegative_intersection = Some(intersection);
            }
        } else {
            if intersection.t > 0. {
                lowest_nonnegative_intersection = Some(intersection);
            }
        }
    }
    lowest_nonnegative_intersection
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene_object::Sphere;

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
        assert_eq!(intersection.object, &sphere);
    }

    #[test]
    fn aggregating_intersections() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(1., &sphere);
        let i2 = Intersection::new(2., &sphere);
        let xs = &[&i1, &i2];
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_all_intersections_have_positive_t() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(1., &sphere);
        let i2 = Intersection::new(2., &sphere);
        let xs = &[&i1, &i2];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i1);
    }

    #[test]
    fn hit_some_intersections_have_negative_t() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(-1., &sphere);
        let i2 = Intersection::new(1., &sphere);
        let xs = &[&i1, &i2];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i2);
    }

    #[test]
    fn hit_all_intersections_have_negative_t() {
        let sphere = Sphere::new();
        let i1 = Intersection::new(-2., &sphere);
        let i2 = Intersection::new(-1., &sphere);
        let xs = &[&i1, &i2];
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
        let xs = &[&i1, &i2, &i3, &i4];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i4);
    }
}
