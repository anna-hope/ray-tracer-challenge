use std::fmt::Debug;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::{
    intersection::{Intersect, Intersection, Ray},
    material::Material,
    Matrix, Result, Tuple,
};

#[derive(Debug, PartialEq)]
pub enum ShapeType {
    Sphere,
    Plane,
    TestShape,
}

pub trait Shape: Intersect + Send + Sync {
    /// Computes the normal vector at the world point.
    fn normal_at(&self, point: Tuple) -> Result<Tuple> {
        let transformation_inverse = self.transformation().inverse()?;
        let local_point = transformation_inverse.clone() * point;
        let local_normal = self.local_normal_at(local_point);
        let mut world_normal = transformation_inverse.transpose() * local_normal;

        // hack to avoid having to find the submatrix of the transformation
        world_normal.w = 0.;
        Ok(world_normal.norm())
    }

    /// Computes the local normal for a given point.
    fn local_normal_at(&self, local_point: Tuple) -> Tuple;

    /// Returns the transformation matrix for this Shape.
    fn transformation(&self) -> Matrix;

    /// Gets object id.
    fn id(&self) -> usize;

    /// Gets object type.
    fn shape_type(&self) -> ShapeType;

    fn material(&self) -> Material {
        Material::default()
    }

    /// Gets an intersection with an arbitrary t for this object.
    /// This is needed primarily for testing, because we can't construct
    /// an intersection with a boxed trait object due to type incompatibility.
    fn arbitrary_intersection(&self, t: f64) -> Intersection;

    /// Sets the object material to the given material.
    /// Needed primarily for testing.
    fn set_material(&mut self, material: Material);
}

impl PartialEq for dyn Shape {
    fn eq(&self, other: &Self) -> bool {
        self.shape_type() == other.shape_type() && self.id() == other.id()
    }
}

impl Debug for dyn Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SceneObject")
            .field("type", &self.shape_type())
            .field("id", &self.id())
            .finish()
    }
}

pub mod sphere {

    use super::*;

    #[derive(Debug, PartialEq, Clone)]
    pub struct Sphere {
        id: usize,
        transformation: Matrix,
        material: Material,
    }

    impl Sphere {
        /// Instantiates a new Sphere with an auto-incrementing id.
        pub fn new() -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            let transformation = Matrix::identity();
            let material = Material::default();
            Self {
                id,
                transformation,
                material,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }

        pub fn with_material(mut self, material: Material) -> Self {
            self.material = material;
            self
        }
    }

    impl Default for Sphere {
        fn default() -> Self {
            Self::new()
        }
    }

    impl Intersect for Sphere {
        /// Calculates the intersection of a sphere and a ray
        /// Returns a Vec of two elements if there is an intersection
        /// (even if it's only in one point, in which case the values would be the same)
        /// or an empty Vec if there is no intersection.
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let ray = ray.transform(self.transformation.inverse()?);

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
                return Ok(vec![]);
            }

            let t1 = (-b - discriminant.sqrt()) / (2. * a);
            let t2 = (-b + discriminant.sqrt()) / (2. * a);

            let i1 = Intersection::new(t1, self);
            let i2 = Intersection::new(t2, self);
            Ok(vec![i1, i2])
        }
    }

    impl Shape for Sphere {
        fn shape_type(&self) -> ShapeType {
            ShapeType::Sphere
        }
        fn id(&self) -> usize {
            self.id
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn local_normal_at(&self, local_point: Tuple) -> Tuple {
            local_point - Tuple::point(0., 0., 0.)
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn arbitrary_intersection(&self, t: f64) -> Intersection {
            Intersection { t, object: self }
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::f64::consts::PI;

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
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 4.0);
            assert_eq!(xs[1].t, 6.0);
        }

        #[test]
        fn ray_intersects_sphere_at_tangent() {
            let ray = Ray::new(Tuple::point(0., 1., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 5.);
            assert_eq!(xs[1].t, 5.);
        }

        #[test]
        fn ray_misses_sphere() {
            let ray = Ray::new(Tuple::point(0., 2., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new();
            let xs = sphere.intersect(&ray).unwrap();
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_originates_inside_sphere() {
            let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, -1.);
            assert_eq!(xs[1].t, 1.);
        }

        #[test]
        fn sphere_is_behind_ray() {
            let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, -6.);
            assert_eq!(xs[1].t, -4.);
        }

        #[test]
        fn intersect_sets_object_on_intersection() {
            let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new();
            let xs = sphere.intersect(&ray).unwrap();
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
            let translation = Matrix::translation(2., 3., 4.);
            let sphere = Sphere::new().with_transformation(translation.clone());
            assert_eq!(sphere.transformation, translation);
        }

        #[test]
        fn intersect_scaled_sphere_with_ray() {
            let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 3.);
            assert_eq!(xs[1].t, 7.);
        }

        #[test]
        fn intersect_translated_sphere_with_ray() {
            let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::new().with_transformation(Matrix::translation(5., 0., 0.));
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 0);
        }

        #[test]
        fn normal_on_sphere_x_axis() {
            let sphere = Sphere::new();
            let normal = sphere.normal_at(Tuple::point(1., 0., 0.)).unwrap();
            assert_eq!(normal, Tuple::vector(1., 0., 0.));
        }

        #[test]
        fn normal_on_sphere_y_axis() {
            let sphere = Sphere::new();
            let normal = sphere.normal_at(Tuple::point(0., 1., 0.)).unwrap();
            assert_eq!(normal, Tuple::vector(0., 1., 0.));
        }

        #[test]
        fn normal_on_sphere_z_axis() {
            let sphere = Sphere::new();
            let normal = sphere.normal_at(Tuple::point(0., 0., 1.)).unwrap();
            assert_eq!(normal, Tuple::vector(0., 0., 1.));
        }

        #[test]
        fn normal_on_sphere_nonaxial() {
            let sphere = Sphere::new();
            let val = 3.0_f64.sqrt() / 3.;
            let normal = sphere.normal_at(Tuple::point(val, val, val)).unwrap();
            assert_eq!(normal, Tuple::vector(val, val, val));
        }

        #[test]
        fn normal_is_normalized_vector() {
            let sphere = Sphere::new();
            let val = 3.0_f64.sqrt() / 3.;
            let normal = sphere.normal_at(Tuple::point(val, val, val)).unwrap();
            assert_eq!(normal, normal.norm());
        }

        #[test]
        fn compute_normal_translated_sphere() {
            let sphere = Sphere::new().with_transformation(Matrix::translation(0., 1., 0.));
            let normal = sphere
                .normal_at(Tuple::point(0., 1.70711, -0.70711))
                .unwrap();
            assert_eq!(normal, Tuple::vector(0., 0.70711, -0.70711));
        }

        #[test]
        fn compute_normal_transformed_sphere() {
            let matrix = Matrix::identity().rotate_z(PI / 5.).scale(1., 0.5, 1.);
            let sphere = Sphere::new().with_transformation(matrix);
            let val = 2.0_f64.sqrt() / 2.;
            let normal = sphere.normal_at(Tuple::point(0., val, -val)).unwrap();
            assert_eq!(normal, Tuple::vector(0., 0.97014, -0.24254));
        }

        #[test]
        fn sphere_has_default_material() {
            let sphere = Sphere::new();
            assert_eq!(sphere.material, Material::default());
        }

        #[test]
        fn sphere_may_be_assigned_material() {
            let mut material = Material::default();
            material.ambient = 1.;
            let sphere = Sphere::new().with_material(material.clone());
            assert_eq!(sphere.material, material);
        }
    }
}

pub mod plane {
    use crate::EPSILON;

    use super::*;

    pub struct Plane {
        id: usize,
        transformation: Matrix,
        material: Material,
    }

    impl Plane {
        pub fn new() -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            let transformation = Matrix::identity();
            let material = Material::default();
            Self {
                id,
                transformation,
                material,
            }
        }

        pub fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            if ray.direction.y.abs() < EPSILON {
                return vec![];
            }

            let t = -ray.origin.y / ray.direction.y;
            vec![Intersection::new(t, self)]
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }

        pub fn with_material(mut self, material: Material) -> Self {
            self.material = material;
            self
        }
    }

    impl Default for Plane {
        fn default() -> Self {
            Self::new()
        }
    }

    impl Intersect for Plane {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let local_ray = ray.transform(self.transformation.inverse()?);
            Ok(self.local_intersect(&local_ray))
        }
    }

    impl Shape for Plane {
        fn arbitrary_intersection(&self, t: f64) -> Intersection {
            Intersection { t, object: self }
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn local_normal_at(&self, _point: Tuple) -> Tuple {
            // Every single point on the plane has the same normal
            Tuple::vector(0., 1., 0.)
        }

        fn id(&self) -> usize {
            self.id
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::Plane
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn normal_of_plane_is_constant_everywhere() {
            let plane = Plane::new();
            let normal_1 = plane.local_normal_at(Tuple::point(0., 0., 0.));
            let normal_2 = plane.local_normal_at(Tuple::point(10., 0., -10.));
            let normal_3 = plane.local_normal_at(Tuple::point(-5., 0., 150.));
            assert_eq!(normal_1, Tuple::vector(0., 1., 0.));
            assert_eq!(normal_2, Tuple::vector(0., 1., 0.));
            assert_eq!(normal_3, Tuple::vector(0., 1., 0.));
        }

        #[test]
        fn intersect_ray_parallel_to_plane() {
            let plane = Plane::new();
            let ray = Ray::new(Tuple::point(0., 10., 0.), Tuple::vector(0., 0., 1.));
            let xs = plane.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn intersect_with_coplanar_ray() {
            let plane = Plane::new();
            let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
            let xs = plane.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_intersecting_plane_from_above() {
            let plane = Plane::new();
            let ray = Ray::new(Tuple::point(0., 1., 0.), Tuple::vector(0., -1., 0.));
            let xs = plane.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0].t, 1.);
            assert_eq!(xs[0].object.shape_type(), ShapeType::Plane);
            assert_eq!(xs[0].object.id(), plane.id);
        }

        #[test]
        fn ray_intersecting_plane_from_below() {
            let plane = Plane::new();
            let ray = Ray::new(Tuple::point(0., -1., 0.), Tuple::vector(0., 1., 0.));
            let xs = plane.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0].t, 1.);
            assert_eq!(xs[0].object.shape_type(), ShapeType::Plane);
            assert_eq!(xs[0].object.id(), plane.id);
        }
    }
}

#[cfg(test)]
mod tests {

    use std::f64::consts::PI;

    use super::*;

    #[derive(Debug)]
    struct TestShape {
        transformation: Matrix,
        material: Material,
    }

    impl TestShape {
        pub fn new() -> Self {
            Self {
                transformation: Matrix::identity(),
                material: Material::default(),
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }

        pub fn with_material(mut self, material: Material) -> Self {
            self.material = material;
            self
        }
    }

    impl Intersect for TestShape {
        fn intersect(&self, _ray: &Ray) -> Result<Vec<Intersection>> {
            unimplemented!()
        }
    }

    impl Shape for TestShape {
        fn id(&self) -> usize {
            0
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::TestShape
        }

        fn arbitrary_intersection(&self, _t: f64) -> Intersection {
            unimplemented!()
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn local_normal_at(&self, local_point: Tuple) -> Tuple {
            Tuple::point(local_point.x, local_point.y, local_point.z)
        }
    }

    #[test]
    fn default_transformation() {
        let shape = TestShape::new();
        assert_eq!(shape.transformation, Matrix::identity());
    }

    #[test]
    fn custom_transformation() {
        let transformation = Matrix::translation(2., 3., 4.);
        let shape = TestShape::new().with_transformation(transformation.clone());
        assert_eq!(shape.transformation, transformation);
    }

    #[test]
    fn default_material() {
        let shape = TestShape::new();
        assert_eq!(shape.material, Material::default());
    }

    #[test]
    fn custom_material() {
        let material = Material {
            ambient: 1.,
            ..Default::default()
        };
        let shape = TestShape::new().with_material(material.clone());
        assert_eq!(shape.material, material);
    }

    #[test]
    fn compute_normal_on_translated_shape() {
        let shape = TestShape::new().with_transformation(Matrix::translation(0., 1., 0.));
        let normal = shape
            .normal_at(Tuple::point(0., 1.70711, -0.70711))
            .unwrap();
        assert_eq!(normal, Tuple::vector(0., 0.70711, -0.70711));
    }

    #[test]
    fn compute_normal_on_transformed_shape() {
        let transformation = Matrix::identity().rotate_z(PI / 2.).scale(1., 0.5, 1.);
        let shape = TestShape::new().with_transformation(transformation);
        let val = 2.0_f64.sqrt() / 2.;
        let normal = shape.normal_at(Tuple::point(0., val, -val)).unwrap();
        assert_eq!(normal, Tuple::vector(0., 0.97014, -0.24254));
    }
}
