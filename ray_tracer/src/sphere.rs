use std::sync::atomic::{AtomicUsize, Ordering};

use crate::{
    intersection::{Intersect, Intersection, Ray},
    material::Material,
    Matrix, Shape, ShapeType, Tuple,
};

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

impl Shape for Sphere {
    fn shape_type(&self) -> ShapeType {
        ShapeType::Sphere
    }
    fn id(&self) -> usize {
        self.id
    }

    fn normal_at(&self, world_point: Tuple) -> Tuple {
        let transformation_inverse = self
            .transformation
            .inverse()
            .expect("Sphere transformation matrix should be invertible");
        let object_point = transformation_inverse.clone() * world_point;
        let object_normal = object_point - Tuple::point(0., 0., 0.);
        let mut world_normal = transformation_inverse.transpose() * object_normal;

        // hack to avoid having to find the submatrix of the transformation
        world_normal.w = 0.;
        world_normal.norm()
    }

    fn material(&self) -> Material {
        self.material
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
        let translation = Matrix::translation(2., 3., 4.);
        let sphere = Sphere::new().with_transformation(translation.clone());
        assert_eq!(sphere.transformation, translation);
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_eq!(xs[0].t, 3.);
        assert_eq!(xs[1].t, 7.);
    }

    #[test]
    fn intersect_translated_sphere_with_ray() {
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let sphere = Sphere::new().with_transformation(Matrix::translation(5., 0., 0.));
        let xs = sphere.intersect(&ray);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn normal_on_sphere_x_axis() {
        let sphere = Sphere::new();
        let normal = sphere.normal_at(Tuple::point(1., 0., 0.));
        assert_eq!(normal, Tuple::vector(1., 0., 0.));
    }

    #[test]
    fn normal_on_sphere_y_axis() {
        let sphere = Sphere::new();
        let normal = sphere.normal_at(Tuple::point(0., 1., 0.));
        assert_eq!(normal, Tuple::vector(0., 1., 0.));
    }

    #[test]
    fn normal_on_sphere_z_axis() {
        let sphere = Sphere::new();
        let normal = sphere.normal_at(Tuple::point(0., 0., 1.));
        assert_eq!(normal, Tuple::vector(0., 0., 1.));
    }

    #[test]
    fn normal_on_sphere_nonaxial() {
        let sphere = Sphere::new();
        let val = 3.0_f64.sqrt() / 3.;
        let normal = sphere.normal_at(Tuple::point(val, val, val));
        assert_eq!(normal, Tuple::vector(val, val, val));
    }

    #[test]
    fn normal_is_normalized_vector() {
        let sphere = Sphere::new();
        let val = 3.0_f64.sqrt() / 3.;
        let normal = sphere.normal_at(Tuple::point(val, val, val));
        assert_eq!(normal, normal.norm());
    }

    #[test]
    fn compute_normal_translated_sphere() {
        let sphere = Sphere::new().with_transformation(Matrix::translation(0., 1., 0.));
        let normal = sphere.normal_at(Tuple::point(0., 1.70711, -0.70711));
        assert_eq!(normal, Tuple::vector(0., 0.70711, -0.70711));
    }

    #[test]
    fn compute_normal_transformed_sphere() {
        let matrix = Matrix::identity().rotate_z(PI / 5.).scale(1., 0.5, 1.);
        let sphere = Sphere::new().with_transformation(matrix);
        let val = 2.0_f64.sqrt() / 2.;
        let normal = sphere.normal_at(Tuple::point(0., val, -val));
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
        let sphere = Sphere::new().with_material(material);
        assert_eq!(sphere.material, material);
    }
}
