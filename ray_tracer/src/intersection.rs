use std::fmt::Debug;

use crate::{
    shape::{Shape, ShapeRef},
    Matrix, Point, Result, Vector, EPSILON,
};

/// Ray.origin is a point, Ray.direction is a vector.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vector,
}

impl Ray {
    pub fn new(origin: Point, direction: Vector) -> Self {
        Self { origin, direction }
    }

    pub fn position(&self, t: f64) -> Point {
        self.origin + self.direction * t
    }

    pub fn transform(&self, matrix: Matrix) -> Self {
        let origin = matrix.clone() * self.origin;
        let direction = matrix * self.direction;
        Self { origin, direction }
    }

    /// Helper method for cylinders and cones.
    ///
    /// Checks to see if the intersection at `t` is within a given radius
    /// from the y axis.
    pub fn check_cap(&self, t: f64, radius: f64) -> bool {
        let x = self.origin.x + t * self.direction.x;
        let z = self.origin.z + t * self.direction.z;

        (x.powi(2) + z.powi(2)) <= radius
    }
}

pub trait Intersect {
    /// Computes the intersections of this ray with the given object.
    fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>>;
}

#[derive(Clone)]
pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a dyn Shape,
    pub point: Point,
    pub eye_vector: Vector,
    pub normal_vector: Vector,
    pub inside: bool,
    pub over_point: Point,
    pub reflect_vector: Vector,
    pub n1: Option<f64>,
    pub n2: Option<f64>,
    pub under_point: Point,
}

impl<'a> Computations<'a> {
    /// Computes the Schlick approximation for Fresnel effect
    /// (i.e. reflection value for transparent surfaces).
    /// Returns 0.0 if either of Computations.n1 and Computations.n2 is None.
    pub fn schlick(&self) -> f64 {
        if let (Some(n1), Some(n2)) = (self.n1, self.n2) {
            // find the cosine of the angle between the eye and the normal vectors
            let mut cos = self.eye_vector.dot(self.normal_vector);

            if n1 > n2 {
                let n = n1 / n2;
                let sin2_t = n.powi(2) * (1. - cos.powi(2));
                if sin2_t > 1. {
                    return 1.;
                }

                // compute cosine of theta_t using trig identity
                let cos_t = (1. - sin2_t).sqrt();

                // when n1 > n2, use cos(theta_t) instead
                cos = cos_t;
            }

            let r0 = ((n1 - n2) / (n1 + n2)).powi(2);
            return r0 + (1. - r0) * (1. - cos).powi(5);
        }

        0.
    }
}

#[derive(Clone, Debug)]
pub struct Intersection {
    pub t: f64,
    pub object: ShapeRef,
    pub u: Option<f64>, // u and v are specific to triangles
    pub v: Option<f64>,
}

impl Intersection {
    pub fn new(t: f64, object: ShapeRef) -> Self {
        Self {
            t,
            object,
            u: None,
            v: None,
        }
    }

    pub fn new_with_uv(t: f64, object: ShapeRef, u: f64, v: f64) -> Self {
        Self {
            t,
            object,
            u: Some(u),
            v: Some(v),
        }
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
        let normal_vector = self.object.normal_at(point, Some(self.clone()))?;

        // if the dot product of normal_vector and eye_vector is negative
        // then they're pointing in (roughly) opposite directions
        let normal_dot_eye = normal_vector.dot(eye_vector);
        let (inside, normal_vector) = if normal_dot_eye < 0. {
            (true, -normal_vector)
        } else {
            (false, normal_vector)
        };

        let reflect_vector = ray.direction.reflect(normal_vector);
        let over_point = point + normal_vector * EPSILON;
        let under_point = point - normal_vector * EPSILON;

        // compute the refractive indices of the materials
        // on either side of a ray-object intersection
        // n1 belongs to the material being exited
        // and n2 belongs to the material being entered
        // (see p. 153 of the book for explanation of the algorithm)
        let mut containers: Vec<&ShapeRef> = vec![];
        let mut n1: Option<f64> = None;
        let mut n2: Option<f64> = None;

        for intersection in intersections {
            if intersection == self {
                if containers.is_empty() {
                    n1 = Some(1.);
                } else if let Some(object) = containers.last() {
                    n1 = Some(object.material().refractive_index);
                }
            }

            let object = &intersection.object;
            if let Some(position) = containers.iter().position(|&x| x == object) {
                containers.remove(position);
            } else {
                containers.push(object);
            }

            if intersection == self {
                if containers.is_empty() {
                    n2 = Some(1.);
                } else if let Some(object) = containers.last() {
                    n2 = Some(object.material().refractive_index);
                }

                break;
            }
        }

        Ok(Computations {
            t: self.t,
            object: &*self.object,
            point,
            eye_vector,
            normal_vector,
            inside,
            over_point,
            reflect_vector,
            n1,
            n2,
            under_point,
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
        equal,
        material::Material,
        shape::{plane::Plane, sphere::Sphere, triangle::Triangle},
        EPSILON,
    };
    use std::sync::Arc;

    #[test]
    fn create_and_query_ray() {
        let origin = Point::new(1., 2., 3.);
        let direction = Vector::new(4., 5., 6.);
        let ray = Ray::new(origin, direction);
        assert_eq!(ray.origin, origin);
        assert_eq!(ray.direction, direction);
    }

    #[test]
    fn compute_point_from_distance() {
        let ray = Ray::new(Point::new(2., 3., 4.), Vector::new(1., 0., 0.));
        assert_eq!(ray.position(0.), Point::new(2., 3., 4.));
        assert_eq!(ray.position(1.), Point::new(3., 3., 4.));
        assert_eq!(ray.position(-1.), Point::new(1., 3., 4.));
        assert_eq!(ray.position(2.5), Point::new(4.5, 3., 4.));
    }

    #[test]
    fn intersection_encapsulates_t_and_object() {
        let sphere = Sphere::default();
        let intersection = Intersection::new(3.5, Arc::new(sphere.clone()));
        assert_eq!(intersection.t, 3.5);
        // assert_eq!(&sphere, intersection.object);
        assert_eq!(intersection.object.id(), sphere.id());
        assert_eq!(intersection.object.shape_type(), sphere.shape_type());
    }

    #[test]
    fn aggregating_intersections() {
        let sphere = Sphere::default();
        let i1 = Intersection::new(1., Arc::new(sphere.clone()));
        let i2 = Intersection::new(2., Arc::new(sphere));
        let xs = &[i1, i2];
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_all_intersections_have_positive_t() {
        let sphere = Sphere::default();
        let i1 = Intersection::new(1., Arc::new(sphere.clone()));
        let i2 = Intersection::new(2., Arc::new(sphere));
        let xs = &[i1.clone(), i2];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i1);
    }

    #[test]
    fn hit_some_intersections_have_negative_t() {
        let sphere = Sphere::default();
        let i1 = Intersection::new(-1., Arc::new(sphere.clone()));
        let i2 = Intersection::new(1., Arc::new(sphere));
        let xs = &[i1, i2.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i2);
    }

    #[test]
    fn hit_all_intersections_have_negative_t() {
        let sphere = Arc::<Sphere>::default();
        let i1 = Intersection::new(-2., sphere.clone());
        let i2 = Intersection::new(-1., sphere);
        let xs = &[i1, i2];
        let intersection = hit(xs);
        assert!(intersection.is_none());
    }

    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let sphere = Arc::<Sphere>::default();
        let i1 = Intersection::new(5., sphere.clone());
        let i2 = Intersection::new(7., sphere.clone());
        let i3 = Intersection::new(-3., sphere.clone());
        let i4 = Intersection::new(2., sphere);
        let xs = &[i1, i2, i3, i4.clone()];
        let intersection = hit(xs).unwrap();
        assert_eq!(intersection, &i4);
    }

    #[test]
    fn translate_ray() {
        let ray = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        let matrix = Matrix::translation(3., 4., 5.);
        let ray2 = ray.transform(matrix);
        assert_eq!(ray2.origin, Point::new(4., 6., 8.));
        assert_eq!(ray2.direction, Vector::new(0., 1., 0.));
    }

    #[test]
    fn scale_ray() {
        let ray = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        let matrix = Matrix::scaling(2., 3., 4.);
        let ray2 = ray.transform(matrix);
        assert_eq!(ray2.origin, Point::new(2., 6., 12.));
        assert_eq!(ray2.direction, Vector::new(0., 3., 0.));
    }

    #[test]
    fn precompute_state_of_intersection() {
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Arc::<Sphere>::default();
        let intersection = Intersection::new(4., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert_eq!(comps.t, intersection.t);
        assert_eq!(comps.object.id(), intersection.object.id());
        assert_eq!(comps.object.shape_type(), intersection.object.shape_type());
        assert_eq!(comps.point, Point::new(0., 0., -1.));
        assert_eq!(comps.eye_vector, Vector::new(0., 0., -1.));
        assert_eq!(comps.normal_vector, Vector::new(0., 0., -1.));
    }

    #[test]
    fn hit_when_intersection_outside() {
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Arc::<Sphere>::default();
        let intersection = Intersection::new(4., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert!(!comps.inside);
    }

    #[test]
    fn hit_when_intersection_inside() {
        let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let shape = Arc::<Sphere>::default();
        let intersection = Intersection::new(1., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert_eq!(comps.point, Point::new(0., 0., 1.));
        assert_eq!(comps.eye_vector, Vector::new(0., 0., -1.));
        assert!(comps.inside);
        assert_eq!(comps.normal_vector, Vector::new(0., 0., -1.));
    }

    #[test]
    fn hit_should_offset_point() {
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(0., 0., 1.)));
        let intersection = Intersection::new(5., shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert!(comps.over_point.z < -EPSILON / 2.);
        assert!(comps.point.z > comps.over_point.z);
    }

    #[test]
    fn precompute_reflection_vector() {
        let shape = Arc::new(Plane::default());
        let val = 2.0_f64 / 2.;
        let ray = Ray::new(Point::new(0., 1., -1.), Vector::new(0., -val, val));
        let intersection = Intersection::new(2.0_f64.sqrt(), shape);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        assert_eq!(comps.reflect_vector, Vector::new(0., val, val));
    }

    #[test]
    fn find_n1_and_n2_at_various_intersections() {
        let a = Arc::new(
            Sphere::glass()
                .with_transformation(Matrix::scaling(2., 2., 2.))
                .with_material(Material {
                    refractive_index: 1.5,
                    ..Default::default()
                }),
        );

        let b = Arc::new(
            Sphere::glass()
                .with_transformation(Matrix::translation(0., 0., -0.25))
                .with_material(Material {
                    refractive_index: 2.,
                    ..Default::default()
                }),
        );

        let c = Arc::new(
            Sphere::glass()
                .with_transformation(Matrix::translation(0., 0., 0.25))
                .with_material(Material {
                    refractive_index: 2.5,
                    ..Default::default()
                }),
        );

        let ray = Ray::new(Point::new(0., 0., -4.), Vector::new(0., 0., 1.));
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

    #[test]
    fn under_point_is_offset_below_surface() {
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Sphere::glass().with_transformation(Matrix::translation(0., 0., 1.));
        let intersection = Intersection::new(5., Arc::new(shape));
        let xs = [intersection.clone()];
        let comps = intersection.prepare_computations(&ray, &xs).unwrap();
        assert!(comps.under_point.z > EPSILON / 2.);
        assert!(comps.point.z < comps.under_point.z);
    }

    #[test]
    fn schlick_approximation_under_total_internal_reflection() {
        let shape = Arc::new(Sphere::glass());
        let val = 2.0_f64 / 2.;
        let ray = Ray::new(Point::new(0., 0., val), Vector::new(0., 1., 0.));
        let xs = &[
            Intersection::new(-val, shape.clone()),
            Intersection::new(val, shape),
        ];
        let comps = xs[1].prepare_computations(&ray, xs).unwrap();
        let reflectance = comps.schlick();
        assert_eq!(reflectance, 1.);
    }

    #[test]
    fn schlick_approximation_with_perpendicular_viewing_angle() {
        let shape = Arc::new(Sphere::glass());
        let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        let xs = &[
            Intersection::new(-1., shape.clone()),
            Intersection::new(1., shape),
        ];
        let comps = xs[1].prepare_computations(&ray, xs).unwrap();
        let reflectance = comps.schlick();
        assert!(equal(reflectance, 0.04));
    }

    #[test]
    fn schlick_approximation_with_small_angle_and_n2_gt_n1() {
        let shape = Arc::new(Sphere::glass());
        let ray = Ray::new(Point::new(0., 0.99, -2.), Vector::new(0., 0., 1.));
        let xs = &[Intersection::new(1.8589, shape)];
        let comps = xs[0].prepare_computations(&ray, xs).unwrap();
        let reflectance = comps.schlick();
        assert!(equal(reflectance, 0.48873));
    }

    #[test]
    fn intersection_with_u_and_v() {
        let shape = Arc::new(Triangle::new(
            Point::new(0., 1., 0.),
            Point::new(-1., 0., 0.),
            Point::new(1., 0., 0.),
            Material::default(),
        ));

        let intersection = Intersection::new_with_uv(3.5, shape, 0.2, 0.4);
        assert_eq!(intersection.u, Some(0.2));
        assert_eq!(intersection.v, Some(0.4));
    }
}
