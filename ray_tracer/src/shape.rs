use std::fmt::Debug;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use lazy_static::lazy_static;
use parking_lot::RwLock;
use slotmap::{DefaultKey, SlotMap};

use crate::{
    intersection::{Intersect, Intersection, Ray},
    material::Material,
    Matrix, Result, Tuple,
};

type ShapeRef = Arc<dyn Shape>;

lazy_static! {
    static ref SHAPES: RwLock<SlotMap<DefaultKey, ShapeRef>> = RwLock::new(SlotMap::new());
}

#[derive(Debug, PartialEq)]
pub enum ShapeType {
    Sphere,
    Plane,
    Cube,
    Cylinder,
    Cone,
    Group,
    TestShape,
}
pub trait ShapeClone {
    fn clone_box(&self) -> Box<dyn Shape>;
}

impl<T> ShapeClone for T
where
    T: 'static + Shape + Clone,
{
    fn clone_box(&self) -> Box<dyn Shape> {
        Box::new(self.clone())
    }
}

pub trait Shape: Intersect + Send + Sync + ShapeClone {
    /// Computes the normal vector at the world point.
    fn normal_at(&self, point: Tuple) -> Result<Tuple> {
        let local_point = self.world_to_object(point)?;
        let local_normal = self.local_normal_at(local_point);
        self.normal_to_world(local_normal)
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

    /// Sets the object material to the given material.
    /// Needed primarily for testing.
    fn set_material(&mut self, material: Material);

    fn parent(&self) -> Option<DefaultKey>;

    fn set_parent(&mut self, parent: DefaultKey);

    fn world_to_object(&self, point: Tuple) -> Result<Tuple> {
        let point = if let Some(parent_key) = self.parent() {
            let shapes = SHAPES.read();
            let parent = Arc::clone(&shapes[parent_key]);
            parent.world_to_object(point)?
        } else {
            point
        };

        Ok(self.transformation().inverse()? * point)
    }

    fn normal_to_world(&self, normal: Tuple) -> Result<Tuple> {
        let mut normal = self.transformation().inverse()?.transpose() * normal;
        normal.w = 0.;
        normal = normal.norm();

        if let Some(parent_key) = self.parent() {
            let shapes = SHAPES.read();
            let parent = Arc::clone(&shapes[parent_key]);
            normal = parent.normal_to_world(normal)?;
        }

        Ok(normal)
    }

    fn add_child(&self, _child: &mut ShapeRef) {
        unimplemented!()
    }
}

impl PartialEq for dyn Shape {
    fn eq(&self, other: &Self) -> bool {
        self.shape_type() == other.shape_type() && self.id() == other.id()
    }
}

impl Debug for dyn Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Shape")
            .field("type", &self.shape_type())
            .field("id", &self.id())
            .finish()
    }
}

impl Clone for Box<dyn Shape> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// Inserts the Arc pointer to the given shape into SHAPES SlotMap.
/// This is typically needed to be called explicitly only for root shapes (`Group`),
/// as child groups/shapes will be inserted into the SlotMap automatically
/// when they are passed to `Group.add_child`.
pub fn insert_shape(shape: Arc<dyn Shape>) -> DefaultKey {
    let mut shapes = SHAPES.write();
    shapes.insert(shape)
}

pub mod sphere {

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Sphere {
        id: usize,
        transformation: Matrix,
        material: Material,
        parent: Option<DefaultKey>,
    }

    impl Sphere {
        /// Instantiates a new Sphere with an auto-incrementing id.
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            Self {
                id,
                transformation,
                material,
                parent,
            }
        }

        /// Instantiates a new glass Sphere
        /// (transparency = 1, refractive_index = 1.5).
        pub fn glass() -> Self {
            let transformation = Matrix::identity();
            let material = Material {
                transparency: 1.,
                refractive_index: 1.5,
                ..Default::default()
            };
            Self {
                transformation,
                material,
                ..Default::default()
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

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
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

            let i1 = Intersection::new(t1, Box::new(self.to_owned()));
            let i2 = Intersection::new(t2, Box::new(self.to_owned()));
            vec![i1, i2]
        }
    }

    impl Default for Sphere {
        fn default() -> Self {
            let transformation = Matrix::identity();
            let material = Material::default();
            Self::new(transformation, material, None)
        }
    }

    impl Intersect for Sphere {
        /// Calculates the intersection of a sphere and a ray
        /// Returns a Vec of two elements if there is an intersection
        /// (even if it's only in one point, in which case the values would be the same)
        /// or an empty Vec if there is no intersection.
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let ray = ray.transform(self.transformation.inverse()?);
            Ok(self.local_intersect(&ray))
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

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::f64::consts::{FRAC_1_SQRT_2, PI};

        #[test]
        fn two_spheres_have_different_ids() {
            let sphere = Sphere::default();
            let sphere2 = Sphere::default();
            assert_ne!(sphere.id, sphere2.id);
        }

        #[test]
        fn ray_intersects_sphere_at_2_points() {
            let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 4.0);
            assert_eq!(xs[1].t, 6.0);
        }

        #[test]
        fn ray_intersects_sphere_at_tangent() {
            let ray = Ray::new(Tuple::point(0., 1., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 5.);
            assert_eq!(xs[1].t, 5.);
        }

        #[test]
        fn ray_misses_sphere() {
            let ray = Ray::new(Tuple::point(0., 2., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_originates_inside_sphere() {
            let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, -1.);
            assert_eq!(xs[1].t, 1.);
        }

        #[test]
        fn sphere_is_behind_ray() {
            let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, -6.);
            assert_eq!(xs[1].t, -4.);
        }

        #[test]
        fn intersect_sets_object_on_intersection() {
            let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].object.id(), sphere.id());
            assert_eq!(xs[1].object.id(), sphere.id());
        }

        #[test]
        fn sphere_default_transformation() {
            let sphere = Sphere::default();
            assert_eq!(sphere.transformation, Matrix::identity());
        }

        #[test]
        fn changing_sphere_transformation() {
            let translation = Matrix::translation(2., 3., 4.);
            let sphere = Sphere::default().with_transformation(translation.clone());
            assert_eq!(sphere.transformation, translation);
        }

        #[test]
        fn intersect_scaled_sphere_with_ray() {
            let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default().with_transformation(Matrix::scaling(2., 2., 2.));
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 3.);
            assert_eq!(xs[1].t, 7.);
        }

        #[test]
        fn intersect_translated_sphere_with_ray() {
            let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
            let sphere = Sphere::default().with_transformation(Matrix::translation(5., 0., 0.));
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 0);
        }

        #[test]
        fn normal_on_sphere_x_axis() {
            let sphere = Sphere::default();
            let normal = sphere.normal_at(Tuple::point(1., 0., 0.)).unwrap();
            assert_eq!(normal, Tuple::vector(1., 0., 0.));
        }

        #[test]
        fn normal_on_sphere_y_axis() {
            let sphere = Sphere::default();
            let normal = sphere.normal_at(Tuple::point(0., 1., 0.)).unwrap();
            assert_eq!(normal, Tuple::vector(0., 1., 0.));
        }

        #[test]
        fn normal_on_sphere_z_axis() {
            let sphere = Sphere::default();
            let normal = sphere.normal_at(Tuple::point(0., 0., 1.)).unwrap();
            assert_eq!(normal, Tuple::vector(0., 0., 1.));
        }

        #[test]
        fn normal_on_sphere_nonaxial() {
            let sphere = Sphere::default();
            let val = 3.0_f64.sqrt() / 3.;
            let normal = sphere.normal_at(Tuple::point(val, val, val)).unwrap();
            assert_eq!(normal, Tuple::vector(val, val, val));
        }

        #[test]
        fn normal_is_normalized_vector() {
            let sphere = Sphere::default();
            let val = 3.0_f64.sqrt() / 3.;
            let normal = sphere.normal_at(Tuple::point(val, val, val)).unwrap();
            assert_eq!(normal, normal.norm());
        }

        #[test]
        fn compute_normal_translated_sphere() {
            let sphere = Sphere::default().with_transformation(Matrix::translation(0., 1., 0.));
            let normal = sphere
                .normal_at(Tuple::point(0., 1.70711, -FRAC_1_SQRT_2))
                .unwrap();
            assert_eq!(normal, Tuple::vector(0., FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
        }

        #[test]
        fn compute_normal_transformed_sphere() {
            let matrix = Matrix::identity().rotate_z(PI / 5.).scale(1., 0.5, 1.);
            let sphere = Sphere::default().with_transformation(matrix);
            let val = 2.0_f64.sqrt() / 2.;
            let normal = sphere.normal_at(Tuple::point(0., val, -val)).unwrap();
            assert_eq!(normal, Tuple::vector(0., 0.97014, -0.24254));
        }

        #[test]
        fn sphere_has_default_material() {
            let sphere = Sphere::default();
            assert_eq!(sphere.material, Material::default());
        }

        #[test]
        fn sphere_may_be_assigned_material() {
            let material = Material {
                ambient: 1.,
                ..Default::default()
            };
            let sphere = Sphere::default().with_material(material.clone());
            assert_eq!(sphere.material, material);
        }

        #[test]
        fn produce_glass_sphere() {
            let sphere = Sphere::glass();
            assert_eq!(sphere.transformation, Matrix::identity());
            assert_eq!(sphere.material.transparency, 1.);
            assert_eq!(sphere.material.refractive_index, 1.5);
        }
    }
}

pub mod plane {
    use crate::EPSILON;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Plane {
        id: usize,
        transformation: Matrix,
        material: Material,
        parent: Option<DefaultKey>,
    }

    impl Plane {
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            Self {
                id,
                transformation,
                material,
                parent,
            }
        }

        pub fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            if ray.direction.y.abs() < EPSILON {
                return vec![];
            }

            let t = -ray.origin.y / ray.direction.y;
            vec![Intersection::new(t, Box::new(self.to_owned()))]
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
            let transformation = Matrix::identity();
            let material = Material::default();
            Self::new(transformation, material, None)
        }
    }

    impl Intersect for Plane {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let local_ray = ray.transform(self.transformation.inverse()?);
            Ok(self.local_intersect(&local_ray))
        }
    }

    impl Shape for Plane {
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

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn normal_of_plane_is_constant_everywhere() {
            let plane = Plane::default();
            let normal_1 = plane.local_normal_at(Tuple::point(0., 0., 0.));
            let normal_2 = plane.local_normal_at(Tuple::point(10., 0., -10.));
            let normal_3 = plane.local_normal_at(Tuple::point(-5., 0., 150.));
            assert_eq!(normal_1, Tuple::vector(0., 1., 0.));
            assert_eq!(normal_2, Tuple::vector(0., 1., 0.));
            assert_eq!(normal_3, Tuple::vector(0., 1., 0.));
        }

        #[test]
        fn intersect_ray_parallel_to_plane() {
            let plane = Plane::default();
            let ray = Ray::new(Tuple::point(0., 10., 0.), Tuple::vector(0., 0., 1.));
            let xs = plane.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn intersect_with_coplanar_ray() {
            let plane = Plane::default();
            let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
            let xs = plane.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_intersecting_plane_from_above() {
            let plane = Plane::default();
            let ray = Ray::new(Tuple::point(0., 1., 0.), Tuple::vector(0., -1., 0.));
            let xs = plane.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0].t, 1.);
            assert_eq!(xs[0].object.shape_type(), ShapeType::Plane);
            assert_eq!(xs[0].object.id(), plane.id);
        }

        #[test]
        fn ray_intersecting_plane_from_below() {
            let plane = Plane::default();
            let ray = Ray::new(Tuple::point(0., -1., 0.), Tuple::vector(0., 1., 0.));
            let xs = plane.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0].t, 1.);
            assert_eq!(xs[0].object.shape_type(), ShapeType::Plane);
            assert_eq!(xs[0].object.id(), plane.id);
        }
    }
}

pub mod cube {
    use std::mem;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Cube {
        id: usize,
        transformation: Matrix,
        material: Material,
        parent: Option<DefaultKey>,
    }

    impl Cube {
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            Self {
                id,
                transformation,
                material,
                parent,
            }
        }

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            // see p. 170 of the book for the intuition behind this algorithm
            let (xt_min, xt_max) = Self::check_axis(ray.origin.x, ray.direction.x);

            // if xt_min == inf or xt_max == -inf, chances are the overall min will be larger than the max
            // which means the ray missed the cube
            if xt_min == f64::INFINITY || xt_max == -f64::INFINITY {
                return vec![];
            }

            let (yt_min, yt_max) = Self::check_axis(ray.origin.y, ray.direction.y);

            // ditto here
            if yt_min == f64::INFINITY || yt_max == -f64::INFINITY {
                return vec![];
            }

            let (zt_min, zt_max) = Self::check_axis(ray.origin.z, ray.direction.z);

            // we know the iterators aren't empty, so ok to unwrap here
            let t_min = *[xt_min, yt_min, zt_min]
                .iter()
                .max_by(|x, y| x.total_cmp(y))
                .unwrap();
            let t_max = *[xt_max, yt_max, zt_max]
                .iter()
                .min_by(|x, y| x.total_cmp(y))
                .unwrap();

            // the minimum t is farther from the ray origin than the maximum t
            // which means the ray misses the cube
            if t_min > t_max {
                return vec![];
            }

            vec![
                Intersection::new(t_min, Box::new(self.clone())),
                Intersection::new(t_max, Box::new(self.clone())),
            ]
        }

        /// Checks to see where the given values intersect the corresponding planes
        /// and returns the minimum and maximum `t` values for each.
        fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
            let t_min_numerator = -1. - origin;
            let t_max_numerator = 1. - origin;

            let mut t_min = t_min_numerator / direction;
            let mut t_max = t_max_numerator / direction;

            if t_min > t_max {
                mem::swap(&mut t_min, &mut t_max);
            }

            (t_min, t_max)
        }
    }

    impl Default for Cube {
        fn default() -> Self {
            Self::new(Matrix::identity(), Material::default(), None)
        }
    }

    impl Shape for Cube {
        fn id(&self) -> usize {
            self.id
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::Cube
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        /// Always picks the component of the point that has the largest absolute value.
        fn local_normal_at(&self, point: Tuple) -> Tuple {
            // we know the array contains elements, so ok to unwrap here
            let max_c = *[point.x.abs(), point.y.abs(), point.z.abs()]
                .iter()
                .max_by(|x, y| x.total_cmp(y))
                .unwrap();

            if max_c == point.x.abs() {
                Tuple::vector(point.x, 0., 0.)
            } else if max_c == point.y.abs() {
                Tuple::vector(0., point.y, 0.)
            } else {
                Tuple::vector(0., 0., point.z)
            }
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material
        }

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }
    }

    impl Intersect for Cube {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let local_ray = ray.transform(self.transformation.inverse()?);
            Ok(self.local_intersect(&local_ray))
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn ray_intersects_cube() {
            let cube = Cube::default();
            let examples = [
                (
                    Tuple::point(5., 0.5, 0.),
                    Tuple::vector(-1., 0., 0.),
                    4.,
                    6.,
                ),
                (
                    Tuple::point(-5., 0.5, 0.),
                    Tuple::vector(1., 0., 0.),
                    4.,
                    6.,
                ),
                (
                    Tuple::point(0.5, 5., 0.),
                    Tuple::vector(0., -1., 0.),
                    4.,
                    6.,
                ),
                (
                    Tuple::point(0.5, -5., 0.),
                    Tuple::vector(0., 1., 0.),
                    4.,
                    6.,
                ),
                (
                    Tuple::point(0.5, 0., 5.),
                    Tuple::vector(0., 0., -1.),
                    4.,
                    6.,
                ),
                (
                    Tuple::point(0.5, 0., -5.),
                    Tuple::vector(0., 0., 1.),
                    4.,
                    6.,
                ),
                (
                    Tuple::point(0., 0.5, 0.),
                    Tuple::vector(0., 0., 1.),
                    -1.,
                    1.,
                ),
            ];

            for (origin, direction, t1, t2) in examples {
                let ray = Ray::new(origin, direction);
                let xs = cube.local_intersect(&ray);
                assert_eq!(xs.len(), 2);
                assert_eq!(xs[0].t, t1);
                assert_eq!(xs[1].t, t2);
            }
        }

        #[test]
        fn ray_misses_cube() {
            let cube = Cube::default();
            let examples = [
                (
                    Tuple::point(-2., 0., 0.),
                    Tuple::vector(0.2673, 0.5345, 0.8018),
                ),
                (
                    Tuple::point(0., -2., 0.),
                    Tuple::vector(0.8018, 0.2673, 0.5345),
                ),
                (
                    Tuple::point(0., 0., -2.),
                    Tuple::vector(0.5345, 0.8018, 0.2673),
                ),
                (Tuple::point(2., 0., 2.), Tuple::vector(0., 0., -1.)),
                (Tuple::point(0., 2., 2.), Tuple::vector(0., -1., 0.)),
                (Tuple::point(2., 2., 0.), Tuple::vector(-1., 0., 0.)),
            ];

            for (origin, direction) in examples {
                let ray = Ray::new(origin, direction);
                let xs = cube.local_intersect(&ray);
                assert_eq!(xs.len(), 0);
            }
        }

        #[test]
        fn normal_on_surface_of_cube() {
            let cube = Cube::default();
            let examples = [
                (Tuple::point(1., 0.5, -0.8), Tuple::vector(1., 0., 0.)),
                (Tuple::point(-1., -0.2, 0.9), Tuple::vector(-1., 0., 0.)),
                (Tuple::point(0.4, 1., -0.1), Tuple::vector(0., 1., 0.)),
                (Tuple::point(0.3, -1., 0.7), Tuple::vector(0., -1., 0.)),
                (Tuple::point(-0.6, 0.3, 1.), Tuple::vector(0., 0., 1.)),
                (Tuple::point(0.4, 0.4, -1.), Tuple::vector(0., 0., -1.)),
                (Tuple::point(1., 1., 1.), Tuple::vector(1., 0., 0.)),
                (Tuple::point(-1., -1., -1.), Tuple::vector(-1., 0., 0.)),
            ];

            for (point, expected_normal) in examples {
                let normal = cube.local_normal_at(point);
                assert_eq!(normal, expected_normal);
            }
        }
    }
}

pub mod cylinder {

    use crate::{equal, EPSILON};
    use std::mem;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Cylinder {
        id: usize,
        transformation: Matrix,
        material: Material,
        minimum: f64,
        maximum: f64,
        closed: bool,
        parent: Option<DefaultKey>,
    }

    impl Cylinder {
        pub fn new(
            transformation: Matrix,
            material: Material,
            minimum: f64,
            maximum: f64,
            closed: bool,
            parent: Option<DefaultKey>,
        ) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            Self {
                id,
                transformation,
                material,
                minimum,
                maximum,
                closed,
                parent,
            }
        }

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            let a = ray.direction.x.powi(2) + ray.direction.z.powi(2);

            // ray is parallel to the y axis
            if equal(a, 0.) {
                return self.intersect_caps(ray);
            }

            let b = 2. * ray.origin.x * ray.direction.x + 2. * ray.origin.z * ray.direction.z;
            let c = ray.origin.x.powi(2) + ray.origin.z.powi(2) - 1.;

            let discriminant = b.powi(2) - 4. * a * c;

            // ray does not intersect the cylinder
            if discriminant < 0. {
                return vec![];
            }

            let mut t0 = (-b - discriminant.sqrt()) / (2. * a);
            let mut t1 = (-b + discriminant.sqrt()) / (2. * a);

            if t0 > t1 {
                mem::swap(&mut t0, &mut t1);
            }

            let mut intersections = vec![];
            let y0 = ray.origin.y + t0 * ray.direction.y;

            if self.minimum < y0 && y0 < self.maximum {
                intersections.push(Intersection::new(t0, Box::new(self.clone())));
            }

            let y1 = ray.origin.y + t1 * ray.direction.y;
            if self.minimum < y1 && y1 < self.maximum {
                intersections.push(Intersection::new(t1, Box::new(self.clone())));
            }

            let mut intersections_caps = self.intersect_caps(ray);
            intersections.append(&mut intersections_caps);

            intersections
        }

        fn intersect_caps(&self, ray: &Ray) -> Vec<Intersection> {
            let mut intersections = vec![];

            // caps only matter if the cylinder is closed, and might possibly be intersected by the ray
            if !self.closed || equal(ray.direction.y, 0.) {
                return intersections;
            }

            // check for an intersection with the lower end cap
            let t = (self.minimum - ray.origin.y) / ray.direction.y;
            if ray.check_cap(t, 1.) {
                intersections.push(Intersection::new(t, Box::new(self.clone())));
            }

            // check for an intersection with the upper end cap
            let t = (self.maximum - ray.origin.y) / ray.direction.y;
            if ray.check_cap(t, 1.) {
                intersections.push(Intersection::new(t, Box::new(self.clone())));
            }

            intersections
        }
    }

    impl Shape for Cylinder {
        fn id(&self) -> usize {
            self.id
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::Cylinder
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn local_normal_at(&self, point: Tuple) -> Tuple {
            // we must check to see which end cap the point corresponds to
            // or see if it lies on the cylinder itself
            // see p.187 of the book for explanation of the algorithm

            // compute the square of the distance from the y axis
            let distance = point.x.powi(2) + point.z.powi(2);

            if distance < 1. && point.y >= self.maximum - EPSILON {
                Tuple::vector(0., 1., 0.)
            } else if distance < 1. && point.y <= self.minimum + EPSILON {
                Tuple::vector(0., -1., 0.)
            } else {
                Tuple::vector(point.x, 0., point.z)
            }
        }

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }
    }

    impl Intersect for Cylinder {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let local_ray = ray.transform(self.transformation.inverse()?);
            Ok(self.local_intersect(&local_ray))
        }
    }

    impl Default for Cylinder {
        fn default() -> Self {
            Self::new(
                Matrix::identity(),
                Material::default(),
                -f64::INFINITY,
                f64::INFINITY,
                false,
                None,
            )
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::equal;

        use super::*;

        #[test]
        fn ray_misses_cylinder() {
            let cylinder = Cylinder::default();
            let examples = [
                (Tuple::point(1., 0., 0.), Tuple::vector(0., 1., 0.)),
                (Tuple::point(0., 0., 0.), Tuple::vector(0., 1., 0.)),
                (Tuple::point(0., 0., -5.), Tuple::vector(1., 1., 1.)),
            ];
            for (origin, direction) in examples {
                let direction = direction.norm();
                let ray = Ray::new(origin, direction);
                let xs = cylinder.local_intersect(&ray);
                assert_eq!(xs.len(), 0);
            }
        }

        #[test]
        fn ray_strikes_cylinder() {
            let cylinder = Cylinder::default();
            let examples = [
                (Tuple::point(1., 0., -5.), Tuple::vector(0., 0., 1.), 5., 5.),
                (Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.), 4., 6.),
                (
                    Tuple::point(0.5, 0., -5.),
                    Tuple::vector(0.1, 1., 1.),
                    6.80798,
                    7.08872,
                ),
            ];

            for (origin, direction, t0, t1) in examples {
                let direction = direction.norm();
                let ray = Ray::new(origin, direction);

                let xs = cylinder.local_intersect(&ray);
                assert_eq!(xs.len(), 2);
                assert!(equal(xs[0].t, t0));
                assert!(equal(xs[1].t, t1));
            }
        }

        #[test]
        fn normal_vector_on_cylinder() {
            let cylinder = Cylinder::default();
            let examples = [
                (Tuple::point(1., 0., 0.), Tuple::vector(1., 0., 0.)),
                (Tuple::point(0., 5., -1.), Tuple::vector(0., 0., -1.)),
                (Tuple::point(0., -2., 1.), Tuple::vector(0., 0., 1.)),
                (Tuple::point(-1., 1., 0.), Tuple::vector(-1., 0., 0.)),
            ];

            for (point, expected_normal) in examples {
                let normal = cylinder.local_normal_at(point);
                assert_eq!(expected_normal, normal);
            }
        }

        #[test]
        fn default_minimum_and_maximum_for_cylinder() {
            let cylinder = Cylinder::default();
            assert_eq!(cylinder.minimum, -f64::INFINITY);
            assert_eq!(cylinder.maximum, f64::INFINITY);
        }

        #[test]
        fn intersecting_constrained_cylinder() {
            let cylinder =
                Cylinder::new(Matrix::identity(), Material::default(), 1., 2., false, None);
            let examples = [
                (Tuple::point(0., 1.5, 0.), Tuple::vector(0.1, 1., 0.), 0),
                (Tuple::point(0., 3., -5.), Tuple::vector(0., 0., 1.), 0),
                (Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.), 0),
                (Tuple::point(0., 2., -5.), Tuple::vector(0., 0., 1.), 0),
                (Tuple::point(0., 1., -5.), Tuple::vector(0., 0., 1.), 0),
                (Tuple::point(0., 1.5, -2.), Tuple::vector(0., 0., 1.), 2),
            ];

            for (point, direction, count) in examples {
                let direction = direction.norm();
                let ray = Ray::new(point, direction);
                let xs = cylinder.local_intersect(&ray);
                assert_eq!(xs.len(), count);
            }
        }

        #[test]
        fn default_closed_cylinder_false() {
            let cylinder = Cylinder::default();
            assert!(!cylinder.closed);
        }

        #[test]
        fn intersecting_caps_of_closed_cylinder() {
            let cylinder =
                Cylinder::new(Matrix::identity(), Material::default(), 1., 2., true, None);
            let examples = [
                (Tuple::point(0., 3., 0.), Tuple::vector(0., -1., 0.), 2),
                (Tuple::point(0., 3., -2.), Tuple::vector(0., -1., 2.), 2),
                (Tuple::point(0., 4., -2.), Tuple::vector(0., -1., 1.), 2),
                (Tuple::point(0., 0., -2.), Tuple::vector(0., 1., 2.), 2),
                (Tuple::point(0., -1., -2.), Tuple::vector(0., 1., 1.), 2),
            ];

            for (point, direction, count) in examples {
                let direction = direction.norm();
                let ray = Ray::new(point, direction);
                let xs = cylinder.local_intersect(&ray);
                assert_eq!(xs.len(), count);
            }
        }

        #[test]
        fn normal_vector_on_cylinders_end_caps() {
            let cylinder =
                Cylinder::new(Matrix::identity(), Material::default(), 1., 2., true, None);
            let examples = [
                (Tuple::point(0., 1., 0.), Tuple::vector(0., -1., 0.)),
                (Tuple::point(0.5, 1., 0.), Tuple::vector(0., -1., 0.)),
                (Tuple::point(0., 1., 0.5), Tuple::vector(0., -1., 0.)),
                (Tuple::point(0., 2., 0.), Tuple::vector(0., 1., 0.)),
                (Tuple::point(0.5, 2., 0.), Tuple::vector(0., 1., 0.)),
                (Tuple::point(0., 2., 0.5), Tuple::vector(0., 1., 0.)),
            ];

            for (point, expected_normal) in examples {
                let normal = cylinder.local_normal_at(point);
                assert_eq!(expected_normal, normal);
            }
        }
    }
}

pub mod cone {
    use std::mem;

    use crate::{equal, EPSILON};

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Cone {
        id: usize,
        transformation: Matrix,
        material: Material,
        minimum: f64,
        maximum: f64,
        closed: bool,
        parent: Option<DefaultKey>,
    }

    impl Cone {
        pub fn new(
            transformation: Matrix,
            material: Material,
            minimum: f64,
            maximum: f64,
            closed: bool,
            parent: Option<DefaultKey>,
        ) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            Self {
                id,
                transformation,
                material,
                minimum,
                maximum,
                closed,
                parent,
            }
        }

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            let a = ray.direction.x.powi(2) - ray.direction.y.powi(2) + ray.direction.z.powi(2);
            let b = 2. * ray.origin.x * ray.direction.x - 2. * ray.origin.y * ray.direction.y
                + 2. * ray.origin.z * ray.direction.z;
            let c = ray.origin.x.powi(2) - ray.origin.y.powi(2) + ray.origin.z.powi(2);

            // ray is parallel to the y axis
            if equal(a, 0.) {
                if self.closed {
                    self.intersect_caps(ray)
                } else {
                    let t = -c / (2. * b);
                    vec![Intersection::new(t, Box::new(self.clone()))]
                }
            } else {
                let discriminant = b.powi(2) - 4. * a * c;

                // ray does not intersect the cylinder
                if discriminant < 0. {
                    return vec![];
                }

                let mut t0 = (-b - discriminant.sqrt()) / (2. * a);
                let mut t1 = (-b + discriminant.sqrt()) / (2. * a);

                if t0 > t1 {
                    mem::swap(&mut t0, &mut t1);
                }

                let mut intersections = vec![];
                let y0 = ray.origin.y + t0 * ray.direction.y;

                if self.minimum < y0 && y0 < self.maximum {
                    intersections.push(Intersection::new(t0, Box::new(self.clone())));
                }

                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection::new(t1, Box::new(self.clone())));
                }

                let mut intersections_caps = self.intersect_caps(ray);
                intersections.append(&mut intersections_caps);

                intersections
            }
        }

        fn intersect_caps(&self, ray: &Ray) -> Vec<Intersection> {
            let mut intersections = vec![];

            // caps only matter if the cone is closed, and might possibly be intersected by the ray
            if !self.closed || equal(ray.direction.y, 0.) {
                return intersections;
            }

            // check for an intersection with the lower end cap
            let t = (self.minimum.abs() - ray.origin.y) / ray.direction.y;
            if ray.check_cap(t, self.minimum.abs()) {
                intersections.push(Intersection::new(t, Box::new(self.clone())));
            }

            // check for an intersection with the upper end cap
            let t = (self.maximum.abs() - ray.origin.y) / ray.direction.y;
            if ray.check_cap(t, self.maximum.abs()) {
                intersections.push(Intersection::new(t, Box::new(self.clone())));
            }

            intersections
        }
    }

    impl Default for Cone {
        fn default() -> Self {
            Self::new(
                Matrix::identity(),
                Material::default(),
                -f64::INFINITY,
                f64::INFINITY,
                false,
                None,
            )
        }
    }

    impl Shape for Cone {
        fn id(&self) -> usize {
            self.id
        }

        fn local_normal_at(&self, point: Tuple) -> Tuple {
            let distance = point.x.powi(2) + point.z.powi(2);

            if distance < 1. && point.y >= self.maximum - EPSILON {
                Tuple::vector(0., 1., 0.)
            } else if distance < 1. && point.y <= self.minimum + EPSILON {
                Tuple::vector(0., -1., 0.)
            } else {
                let y = (point.x.powi(2) + point.z.powi(2)).sqrt();
                let y = if point.y > 0. { -y } else { y };
                Tuple::vector(point.x, y, point.z)
            }
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::Cone
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }
    }

    impl Intersect for Cone {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let local_ray = ray.transform(self.transformation.inverse()?);
            Ok(self.local_intersect(&local_ray))
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn intersecting_cone_with_ray() {
            let shape = Cone::default();
            let examples = [
                (Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.), 5., 5.),
                (
                    Tuple::point(0., 0., -5.),
                    Tuple::vector(1., 1., 1.),
                    8.66025,
                    8.66025,
                ),
                (
                    Tuple::point(1., 1., -5.),
                    Tuple::vector(-0.5, -1., 1.),
                    4.55006,
                    49.44994,
                ),
            ];

            for (origin, direction, t0, t1) in examples {
                let direction = direction.norm();
                let ray = Ray::new(origin, direction);
                let xs = shape.local_intersect(&ray);

                assert_eq!(xs.len(), 2);
                assert!(equal(xs[0].t, t0));
                assert!(equal(xs[1].t, t1));
            }
        }

        #[test]
        fn intersect_cone_with_ray_parallel_to_one_of_its_halves() {
            let shape = Cone::default();
            let direction = Tuple::vector(0., 1., 1.).norm();
            let ray = Ray::new(Tuple::point(0., 0., -1.), direction);
            let xs = shape.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert!(equal(xs[0].t, 0.35355));
        }

        #[test]
        fn intersecting_cones_end_caps() {
            let shape = Cone::new(
                Matrix::identity(),
                Material::default(),
                -0.5,
                0.5,
                true,
                None,
            );
            let examples = [
                (Tuple::point(0., 0., -5.), Tuple::vector(0., 1., 0.), 0),
                (Tuple::point(0., 0., -0.25), Tuple::vector(0., 1., 1.), 2),
                (Tuple::point(0., 0., -0.25), Tuple::vector(0., 1., 0.), 4),
            ];

            for (origin, direction, count) in examples {
                let direction = direction.norm();
                let ray = Ray::new(origin, direction);
                let xs = shape.local_intersect(&ray);
                assert_eq!(xs.len(), count);
            }
        }

        #[test]
        fn compute_normal_vector_on_cone() {
            let shape = Cone::default();
            let examples = [
                (Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 0.)),
                (
                    Tuple::point(1., 1., 1.),
                    Tuple::vector(1., -2.0_f64.sqrt(), 1.),
                ),
                (Tuple::point(-1., -1., 0.), Tuple::vector(-1., 1., 0.)),
            ];

            for (point, expected_normal) in examples {
                let normal = shape.local_normal_at(point);
                assert_eq!(expected_normal, normal);
            }
        }
    }
}

pub mod group {
    type GroupChildren = Arc<RwLock<Vec<DefaultKey>>>;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Group {
        id: usize,
        transformation: Matrix,
        children: GroupChildren,
        parent: Option<DefaultKey>,
    }

    impl Group {
        pub fn new(
            transformation: Matrix,
            children: GroupChildren,
            parent: Option<DefaultKey>,
        ) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);

            Self {
                id,
                transformation,
                children,
                parent,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }

        fn local_intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let children = self.children.read();
            let mut intersections = vec![];
            let shapes = SHAPES.read();
            for child_key in children.iter() {
                let child = shapes.get(*child_key).unwrap();
                let mut child_intersections = child.intersect(ray)?;
                intersections.append(&mut child_intersections);
            }

            intersections.sort_unstable_by(|a, b| a.t.total_cmp(&b.t));
            Ok(intersections)
        }
    }

    impl Default for Group {
        fn default() -> Self {
            let children = Arc::new(RwLock::new(vec![]));
            Self::new(Matrix::identity(), children, None)
        }
    }

    impl Shape for Group {
        fn id(&self) -> usize {
            self.id
        }

        fn local_normal_at(&self, _: Tuple) -> Tuple {
            unimplemented!()
        }

        fn material(&self) -> Material {
            unimplemented!()
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn set_material(&mut self, _: Material) {
            unimplemented!()
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::Group
        }

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }

        fn add_child(&self, child: &mut ShapeRef) {
            {
                let shapes = SHAPES.read();
                let mut self_key: Option<DefaultKey> = None;
                for (key, value) in shapes.iter() {
                    if value.id() == self.id() && value.shape_type() == ShapeType::Group {
                        self_key = Some(key);
                    }
                }

                let self_key = self_key.expect("The group must be in SHAPES");
                let child_shape =
                    Arc::get_mut(child).expect("Must have a unique reference to the child");
                child_shape.set_parent(self_key);
            }

            let child_key = insert_shape(Arc::clone(child));
            let mut children = self.children.write();
            children.push(child_key);
        }
    }

    impl Intersect for Group {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let ray = ray.transform(self.transformation.inverse()?);
            self.local_intersect(&ray)
        }
    }

    impl Drop for Group {
        fn drop(&mut self) {
            println!("Dropping group {}", self.id);
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::shape::sphere::Sphere;

        use super::*;

        #[test]
        fn create_new_group() {
            let group = Group::default();
            assert_eq!(group.transformation, Matrix::identity());
            let children = group.children.read();
            assert!(children.is_empty());
        }

        #[test]
        fn add_child_to_group() {
            // a somewhat convoluted test to make sure that we can add shapes to a group
            // that the group then has those shapes as their children
            // that the parent(s) of the shapes refer to the same group
            // and that the parent(s) of those shapes (being the same group)
            // contain those shapes
            let shape = super::super::tests::TestShape::default();

            let group = Arc::new(Group::default());
            let group_clone: ShapeRef = group.clone() as ShapeRef;
            insert_shape(group_clone);

            let mut shape: ShapeRef = Arc::new(shape);
            group.add_child(&mut shape);

            let children = group.children.read();
            assert!(!children.is_empty());

            let shapes = SHAPES.read();
            let child_key = *children.first().unwrap();
            let child = shapes.get(child_key).unwrap();

            assert_eq!(shape.id(), child.id());
            assert_eq!(shape.shape_type(), child.shape_type());

            let parent_key = child.parent().unwrap();
            let parent = shapes.get(parent_key).unwrap();
            assert_eq!(parent.id(), group.id);
        }

        #[test]
        fn intersect_ray_with_empty_group() {
            let group = Group::default();
            let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));
            let xs = group.local_intersect(&ray).unwrap();
            assert!(xs.is_empty());
        }

        #[test]
        fn intersecting_ray_with_nonempty_group() {
            let s1 = Sphere::default();
            let s2 = Sphere::default().with_transformation(Matrix::translation(0., 0., -3.));
            let s3 = Sphere::default().with_transformation(Matrix::translation(5., 0., 0.));

            let group = Arc::new(Group::default());
            let group_clone: ShapeRef = Arc::clone(&group) as ShapeRef;
            insert_shape(group_clone);

            let mut s1: ShapeRef = Arc::new(s1);
            let mut s2: ShapeRef = Arc::new(s2);
            let mut s3: ShapeRef = Arc::new(s3);

            group.add_child(&mut s1);
            group.add_child(&mut s2);
            group.add_child(&mut s3);

            let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
            let xs = group.local_intersect(&ray).unwrap();

            assert_eq!(xs.len(), 4);
            assert_eq!(xs[0].object.id(), s2.id());
            assert_eq!(xs[1].object.id(), s2.id());
            assert_eq!(xs[2].object.id(), s1.id());
            assert_eq!(xs[3].object.id(), s1.id());
        }

        #[test]
        fn intersect_transformed_group() {
            let shape = Sphere::default().with_transformation(Matrix::translation(5., 0., 0.));

            let group = Arc::new(Group::default().with_transformation(Matrix::scaling(2., 2., 2.)));
            let group_clone = Arc::clone(&group) as ShapeRef;
            insert_shape(group_clone);

            let mut shape: ShapeRef = Arc::new(shape);
            group.add_child(&mut shape);

            let ray = Ray::new(Tuple::point(10., 0., -10.), Tuple::vector(0., 0., 1.));
            let xs = group.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
        }
    }
}

#[cfg(test)]
mod tests {

    use std::f64::consts::{FRAC_1_SQRT_2, PI};

    use super::{group::Group, sphere::Sphere, *};

    #[derive(Debug, Clone)]
    pub struct TestShape {
        transformation: Matrix,
        material: Material,
        parent: Option<DefaultKey>,
    }

    impl TestShape {
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            Self {
                transformation,
                material,
                parent,
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

    impl Default for TestShape {
        fn default() -> Self {
            Self::new(Matrix::identity(), Material::default(), None)
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

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn local_normal_at(&self, local_point: Tuple) -> Tuple {
            Tuple::point(local_point.x, local_point.y, local_point.z)
        }

        fn parent(&self) -> Option<DefaultKey> {
            self.parent
        }

        fn set_parent(&mut self, parent: DefaultKey) {
            self.parent = Some(parent);
        }
    }

    #[test]
    fn default_transformation() {
        let shape = TestShape::default();
        assert_eq!(shape.transformation, Matrix::identity());
    }

    #[test]
    fn custom_transformation() {
        let transformation = Matrix::translation(2., 3., 4.);
        let shape = TestShape::default().with_transformation(transformation.clone());
        assert_eq!(shape.transformation, transformation);
    }

    #[test]
    fn default_material() {
        let shape = TestShape::default();
        assert_eq!(shape.material, Material::default());
    }

    #[test]
    fn custom_material() {
        let material = Material {
            ambient: 1.,
            ..Default::default()
        };
        let shape = TestShape::default().with_material(material.clone());
        assert_eq!(shape.material, material);
    }

    #[test]
    fn compute_normal_on_translated_shape() {
        let shape = TestShape::default().with_transformation(Matrix::translation(0., 1., 0.));
        let normal = shape
            .normal_at(Tuple::point(0., 1.70711, -FRAC_1_SQRT_2))
            .unwrap();
        assert_eq!(normal, Tuple::vector(0., FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn compute_normal_on_transformed_shape() {
        let transformation = Matrix::identity().rotate_z(PI / 2.).scale(1., 0.5, 1.);
        let shape = TestShape::default().with_transformation(transformation);
        let val = 2.0_f64.sqrt() / 2.;
        let normal = shape.normal_at(Tuple::point(0., val, -val)).unwrap();
        assert_eq!(normal, Tuple::vector(0., 0.97014, -0.24254));
    }

    #[test]
    fn shape_has_parent_field() {
        let shape = TestShape::default();
        assert!(shape.parent.is_none());
    }

    #[test]
    fn convert_point_from_world_to_object_space() {
        let mut sphere: ShapeRef =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(5., 0., 0.)));

        let group1 = Arc::new(Group::default().with_transformation(Matrix::rotation_y(PI / 2.)));
        let mut group2: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::scaling(2., 2., 2.)));

        let group1_clone: ShapeRef = Arc::clone(&group1) as ShapeRef;
        insert_shape(group1_clone);

        group1.add_child(&mut group2);
        group2.add_child(&mut sphere);

        let point = sphere.world_to_object(Tuple::point(-2., 0., -10.)).unwrap();
        assert_eq!(point, Tuple::point(0., 0., -1.));
    }

    #[test]
    fn convert_normal_from_object_to_world_space() {
        let group1 = Arc::new(Group::default().with_transformation(Matrix::rotation_y(PI / 2.)));
        let group1_clone = Arc::clone(&group1) as ShapeRef;
        insert_shape(group1_clone);

        let mut group2: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::scaling(1., 2., 3.)));
        group1.add_child(&mut group2);

        let mut sphere: ShapeRef =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(5., 0., 0.)));
        group2.add_child(&mut sphere);

        let val = 3.0_f64.sqrt() / 3.;
        let normal = sphere
            .normal_to_world(Tuple::vector(val, val, val))
            .unwrap();

        // need more precise values (5 significant figures)
        // than the book provides
        // cf. book values (0.2857, 0.4286, -0.8571)
        assert_eq!(normal, Tuple::vector(0.28571, 0.42857, -0.85714));
    }

    #[test]
    fn find_normal_on_child_object() {
        let group1: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::rotation_y(PI / 2.)));
        insert_shape(Arc::clone(&group1));

        let mut group2: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::scaling(1., 2., 3.)));
        group1.add_child(&mut group2);

        let mut sphere: ShapeRef =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(5., 0., 0.)));
        group2.add_child(&mut sphere);

        let normal = sphere
            .normal_at(Tuple::point(1.7321, 1.1547, -5.5774))
            .unwrap();
        // cf. book values (0.2857, 0.4286, -0.8571)
        assert_eq!(normal, Tuple::vector(0.2857, 0.42854, -0.85716));
    }
}
