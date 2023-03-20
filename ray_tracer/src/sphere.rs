use std::sync::atomic::{AtomicUsize, Ordering};

use crate::{
    intersection::{Intersect, Intersection, Ray},
    Id, Matrix, ObjectType, SceneObject, SceneObjectType, Tuple,
};

#[derive(Debug, PartialEq, Clone)]
pub struct Sphere {
    id: usize,
    pub transformation: Matrix,
}

impl Sphere {
    /// Instantiates a new Sphere with an auto-incrementing id.
    pub fn new() -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(1);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let transformation = Matrix::identity();
        Self { id, transformation }
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self::new()
    }
}

impl Id for Sphere {
    fn id(&self) -> usize {
        self.id
    }
}

impl ObjectType for Sphere {
    fn object_type(&self) -> SceneObjectType {
        SceneObjectType::Sphere
    }
}

impl Intersect for Sphere {
    /// Calculates the intersection of a sphere and a ray
    /// Returns a Vec of two elements if there is an intersection
    /// (even if it's only in one point, in which case the values would be the same)
    /// or an empty Vec if there is no intersection.
    fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let ray = ray.transform(
            self.transformation
                .inverse()
                .expect("Sphere transformation should be invertible"),
        );

        // the vector from the sphere's center, to the ray origin
        // (the sphere is centered at the world origin)
        // (subtracting a point from a point gives us a vector)
        let sphere_to_ray = ray.origin - Tuple::point(0., 0., 0.);

        // ok to unwrap() here, since we know for sure sphere_to_ray is a vector
        let a = ray.direction.dot(&ray.direction).unwrap();
        let b = 2. * ray.direction.dot(&sphere_to_ray).unwrap();
        let c = sphere_to_ray.dot(&sphere_to_ray).unwrap() - 1.;

        let discriminant = b.powi(2) - 4. * a * c;
        if discriminant < 0. {
            return vec![];
        }

        let t1 = (-b - discriminant.sqrt()) / (2. * a);
        let t2 = (-b + discriminant.sqrt()) / (2. * a);

        let i1 = Intersection::new(t1, self);
        let i2 = Intersection::new(t2, self);
        vec![i1, i2]
    }
}

impl SceneObject for Sphere {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_spheres_have_different_ids() {
        let sphere = Sphere::new();
        let sphere2 = Sphere::new();
        assert_ne!(sphere.id, sphere2.id);
    }

    #[test]
    fn ray_intersects_sphere_at_2_points() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new();
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 4.0);
        assert_eq!(xs[1].t, 6.0);
    }

    #[test]
    fn ray_intersects_sphere_at_tangent() {
        let ray = Ray::new(Tuple::point(0., 1., -5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new();
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 5.);
        assert_eq!(xs[1].t, 5.);
    }

    #[test]
    fn ray_misses_sphere() {
        let ray = Ray::new(Tuple::point(0., 2., -5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new();
        let xs = sphere.intersect(&ray);
        assert!(xs.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new();
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, -1.);
        assert_eq!(xs[1].t, 1.);
    }

    #[test]
    fn sphere_is_behind_ray() {
        let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new();
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, -6.);
        assert_eq!(xs[1].t, -4.);
    }

    #[test]
    fn intersect_sets_object_on_intersection() {
        let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new();
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].object.id(), sphere.id());
        assert_eq!(xs[1].object.id(), sphere.id());
    }

    #[test]
    fn sphere_default_transformation() {
        let sphere = Sphere::new();
        assert_eq!(sphere.transformation, Matrix::identity());
    }

    #[test]
    fn changing_sphere_transformation() {
        let mut sphere = Sphere::new();
        let translation = Matrix::translation(2., 3., 4.);
        sphere.transformation = translation.clone();
        assert_eq!(sphere.transformation, translation);
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut sphere = Sphere::new();
        sphere.transformation = Matrix::scaling(2., 2., 2.);
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 3.);
        assert_eq!(xs[1].t, 7.);
    }

    #[test]
    fn intersect_translated_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let mut sphere = Sphere::new();
        sphere.transformation = Matrix::translation(5., 0., 0.);
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 0);
    }
}
