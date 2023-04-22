use std::any::Any;
use std::fmt::Debug;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use lazy_static::lazy_static;
use parking_lot::RwLock;
use slotmap::{DefaultKey, SlotMap};

use crate::{
    intersection::{Intersect, Intersection, Ray},
    material::Material,
    Matrix, Point, Result, Vector,
};

use self::group::Group;

pub type ShapeRef = Arc<dyn Shape>;
type MaybeKeyRef = Arc<RwLock<Option<DefaultKey>>>;

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
    Triangle,
    Group,
    TestShape,
}

pub trait Shape: Intersect + Send + Sync {
    /// Computes the normal vector at the world point.
    /// Typically does not need to be implemented for concrete types.
    fn normal_at(&self, point: Point) -> Result<Vector> {
        let local_point = self.world_to_object(point)?;
        let local_normal = self.local_normal_at(local_point);
        self.normal_to_world(local_normal)
    }

    /// Computes the local normal for a given point.
    fn local_normal_at(&self, local_point: Point) -> Vector;

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

    fn parent_key(&self) -> Option<DefaultKey>;

    fn parent(&self) -> Option<ShapeRef> {
        let parent_key = self.parent_key()?;
        let shapes = SHAPES.read();
        let parent = shapes.get(parent_key)?;
        Some(Arc::clone(parent))
    }

    fn set_parent(&self, parent: DefaultKey);

    /// Converts the given point from world space to object space.
    /// Typically does not need to be implemented for concrete types.
    fn world_to_object(&self, point: Point) -> Result<Point> {
        let point = if let Some(parent) = self.parent() {
            parent.world_to_object(point)?
        } else {
            point
        };

        Ok(self.transformation().inverse()? * point)
    }

    /// Converts the given normal to world space.
    /// Typically does not need to be implemented for concrete types.
    fn normal_to_world(&self, normal: Vector) -> Result<Vector> {
        let mut normal = self.transformation().inverse()?.transpose() * normal;
        normal = normal.norm();

        if let Some(parent) = self.parent() {
            normal = parent.normal_to_world(normal)?;
        }

        Ok(normal)
    }

    /// Upcast self to Any. This is needed primarily for testing,
    /// since in some cases we need to get back to original types.
    fn as_any(&self) -> &dyn Any;

    /// Downcasts this shape to a Group, if it is actually a Group.
    /// Returns None otherwise. Should not typically be implemented for concrete types.
    fn as_group(&self) -> Option<&Group> {
        self.as_any().downcast_ref::<Group>()
    }

    /// Gets the key of this shape from the SHAPES store, if the shape is in it.
    /// This is an O(n) operation, where n is the size of the SHAPES store,
    /// since we need to iterate over all the shapes in the store
    /// and compare them with this shape for equality.
    /// Should not typically be implemented for concrete types.
    fn key(&self) -> Option<DefaultKey> {
        let shapes = SHAPES.read();
        for (key, shape) in shapes.iter() {
            if shape.id() == self.id() && shape.shape_type() == self.shape_type() {
                return Some(key);
            }
        }
        None
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

/// Inserts the Arc pointer to the given shape into SHAPES SlotMap.
/// This is typically needed to be called explicitly only for root shapes (`Group`),
/// as child groups/shapes will be inserted into the SlotMap automatically
/// when they are passed to `Group.add_child`.
pub fn register_shape(shape: ShapeRef) -> DefaultKey {
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
        parent: MaybeKeyRef,
    }

    impl Sphere {
        /// Instantiates a new Sphere with an auto-incrementing id.
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);

            let parent = Arc::new(RwLock::new(parent));
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
            let sphere_to_ray = ray.origin - Point::new(0., 0., 0.);

            let a = ray.direction.dot(ray.direction);
            let b = 2. * ray.direction.dot(sphere_to_ray);
            let c = sphere_to_ray.dot(sphere_to_ray) - 1.;

            let discriminant = b.powi(2) - 4. * a * c;
            if discriminant < 0. {
                return vec![];
            }

            let t1 = (-b - discriminant.sqrt()) / (2. * a);
            let t2 = (-b + discriminant.sqrt()) / (2. * a);

            let i1 = Intersection::new(t1, Arc::new(self.to_owned()));
            let i2 = Intersection::new(t2, Arc::new(self.to_owned()));
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

        fn local_normal_at(&self, local_point: Point) -> Vector {
            local_point - Point::new(0., 0., 0.)
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
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
            let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 4.0);
            assert_eq!(xs[1].t, 6.0);
        }

        #[test]
        fn ray_intersects_sphere_at_tangent() {
            let ray = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 5.);
            assert_eq!(xs[1].t, 5.);
        }

        #[test]
        fn ray_misses_sphere() {
            let ray = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_originates_inside_sphere() {
            let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, -1.);
            assert_eq!(xs[1].t, 1.);
        }

        #[test]
        fn sphere_is_behind_ray() {
            let ray = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default();
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, -6.);
            assert_eq!(xs[1].t, -4.);
        }

        #[test]
        fn intersect_sets_object_on_intersection() {
            let ray = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
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
            let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default().with_transformation(Matrix::scaling(2., 2., 2.));
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
            assert_eq!(xs[0].t, 3.);
            assert_eq!(xs[1].t, 7.);
        }

        #[test]
        fn intersect_translated_sphere_with_ray() {
            let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
            let sphere = Sphere::default().with_transformation(Matrix::translation(5., 0., 0.));
            let xs = sphere.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 0);
        }

        #[test]
        fn normal_on_sphere_x_axis() {
            let sphere = Sphere::default();
            let normal = sphere.normal_at(Point::new(1., 0., 0.)).unwrap();
            assert_eq!(normal, Vector::new(1., 0., 0.));
        }

        #[test]
        fn normal_on_sphere_y_axis() {
            let sphere = Sphere::default();
            let normal = sphere.normal_at(Point::new(0., 1., 0.)).unwrap();
            assert_eq!(normal, Vector::new(0., 1., 0.));
        }

        #[test]
        fn normal_on_sphere_z_axis() {
            let sphere = Sphere::default();
            let normal = sphere.normal_at(Point::new(0., 0., 1.)).unwrap();
            assert_eq!(normal, Vector::new(0., 0., 1.));
        }

        #[test]
        fn normal_on_sphere_nonaxial() {
            let sphere = Sphere::default();
            let val = 3.0_f64.sqrt() / 3.;
            let normal = sphere.normal_at(Point::new(val, val, val)).unwrap();
            assert_eq!(normal, Vector::new(val, val, val));
        }

        #[test]
        fn normal_is_normalized_vector() {
            let sphere = Sphere::default();
            let val = 3.0_f64.sqrt() / 3.;
            let normal = sphere.normal_at(Point::new(val, val, val)).unwrap();
            assert_eq!(normal, normal.norm());
        }

        #[test]
        fn compute_normal_translated_sphere() {
            let sphere = Sphere::default().with_transformation(Matrix::translation(0., 1., 0.));
            let normal = sphere
                .normal_at(Point::new(0., 1.70711, -FRAC_1_SQRT_2))
                .unwrap();
            assert_eq!(normal, Vector::new(0., FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
        }

        #[test]
        fn compute_normal_transformed_sphere() {
            let matrix = Matrix::identity().rotate_z(PI / 5.).scale(1., 0.5, 1.);
            let sphere = Sphere::default().with_transformation(matrix);
            let val = 2.0_f64.sqrt() / 2.;
            let normal = sphere.normal_at(Point::new(0., val, -val)).unwrap();
            assert_eq!(normal, Vector::new(0., 0.97014, -0.24254));
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
        parent: MaybeKeyRef,
    }

    impl Plane {
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            let parent = Arc::new(RwLock::new(parent));
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
            vec![Intersection::new(t, Arc::new(self.to_owned()))]
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

        fn local_normal_at(&self, _point: Point) -> Vector {
            // Every single point on the plane has the same normal
            Vector::new(0., 1., 0.)
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

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn normal_of_plane_is_constant_everywhere() {
            let plane = Plane::default();
            let normal_1 = plane.local_normal_at(Point::new(0., 0., 0.));
            let normal_2 = plane.local_normal_at(Point::new(10., 0., -10.));
            let normal_3 = plane.local_normal_at(Point::new(-5., 0., 150.));
            assert_eq!(normal_1, Vector::new(0., 1., 0.));
            assert_eq!(normal_2, Vector::new(0., 1., 0.));
            assert_eq!(normal_3, Vector::new(0., 1., 0.));
        }

        #[test]
        fn intersect_ray_parallel_to_plane() {
            let plane = Plane::default();
            let ray = Ray::new(Point::new(0., 10., 0.), Vector::new(0., 0., 1.));
            let xs = plane.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn intersect_with_coplanar_ray() {
            let plane = Plane::default();
            let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
            let xs = plane.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_intersecting_plane_from_above() {
            let plane = Plane::default();
            let ray = Ray::new(Point::new(0., 1., 0.), Vector::new(0., -1., 0.));
            let xs = plane.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0].t, 1.);
            assert_eq!(xs[0].object.shape_type(), ShapeType::Plane);
            assert_eq!(xs[0].object.id(), plane.id);
        }

        #[test]
        fn ray_intersecting_plane_from_below() {
            let plane = Plane::default();
            let ray = Ray::new(Point::new(0., -1., 0.), Vector::new(0., 1., 0.));
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
        parent: MaybeKeyRef,
    }

    impl Cube {
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            let parent = Arc::new(RwLock::new(parent));
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
                Intersection::new(t_min, Arc::new(self.clone())),
                Intersection::new(t_max, Arc::new(self.clone())),
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
        fn local_normal_at(&self, point: Point) -> Vector {
            // we know the array contains elements, so ok to unwrap here
            let max_c = *[point.x.abs(), point.y.abs(), point.z.abs()]
                .iter()
                .max_by(|x, y| x.total_cmp(y))
                .unwrap();

            if max_c == point.x.abs() {
                Vector::new(point.x, 0., 0.)
            } else if max_c == point.y.abs() {
                Vector::new(0., point.y, 0.)
            } else {
                Vector::new(0., 0., point.z)
            }
        }

        fn material(&self) -> Material {
            self.material.clone()
        }

        fn set_material(&mut self, material: Material) {
            self.material = material
        }

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
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
                (Point::new(5., 0.5, 0.), Vector::new(-1., 0., 0.), 4., 6.),
                (Point::new(-5., 0.5, 0.), Vector::new(1., 0., 0.), 4., 6.),
                (Point::new(0.5, 5., 0.), Vector::new(0., -1., 0.), 4., 6.),
                (Point::new(0.5, -5., 0.), Vector::new(0., 1., 0.), 4., 6.),
                (Point::new(0.5, 0., 5.), Vector::new(0., 0., -1.), 4., 6.),
                (Point::new(0.5, 0., -5.), Vector::new(0., 0., 1.), 4., 6.),
                (Point::new(0., 0.5, 0.), Vector::new(0., 0., 1.), -1., 1.),
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
                (Point::new(-2., 0., 0.), Vector::new(0.2673, 0.5345, 0.8018)),
                (Point::new(0., -2., 0.), Vector::new(0.8018, 0.2673, 0.5345)),
                (Point::new(0., 0., -2.), Vector::new(0.5345, 0.8018, 0.2673)),
                (Point::new(2., 0., 2.), Vector::new(0., 0., -1.)),
                (Point::new(0., 2., 2.), Vector::new(0., -1., 0.)),
                (Point::new(2., 2., 0.), Vector::new(-1., 0., 0.)),
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
                (Point::new(1., 0.5, -0.8), Vector::new(1., 0., 0.)),
                (Point::new(-1., -0.2, 0.9), Vector::new(-1., 0., 0.)),
                (Point::new(0.4, 1., -0.1), Vector::new(0., 1., 0.)),
                (Point::new(0.3, -1., 0.7), Vector::new(0., -1., 0.)),
                (Point::new(-0.6, 0.3, 1.), Vector::new(0., 0., 1.)),
                (Point::new(0.4, 0.4, -1.), Vector::new(0., 0., -1.)),
                (Point::new(1., 1., 1.), Vector::new(1., 0., 0.)),
                (Point::new(-1., -1., -1.), Vector::new(-1., 0., 0.)),
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
        parent: MaybeKeyRef,
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
            let parent = Arc::new(RwLock::new(parent));
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
                intersections.push(Intersection::new(t0, Arc::new(self.clone())));
            }

            let y1 = ray.origin.y + t1 * ray.direction.y;
            if self.minimum < y1 && y1 < self.maximum {
                intersections.push(Intersection::new(t1, Arc::new(self.clone())));
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
                intersections.push(Intersection::new(t, Arc::new(self.clone())));
            }

            // check for an intersection with the upper end cap
            let t = (self.maximum - ray.origin.y) / ray.direction.y;
            if ray.check_cap(t, 1.) {
                intersections.push(Intersection::new(t, Arc::new(self.clone())));
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

        fn local_normal_at(&self, point: Point) -> Vector {
            // we must check to see which end cap the point corresponds to
            // or see if it lies on the cylinder itself
            // see p.187 of the book for explanation of the algorithm

            // compute the square of the distance from the y axis
            let distance = point.x.powi(2) + point.z.powi(2);

            if distance < 1. && point.y >= self.maximum - EPSILON {
                Vector::new(0., 1., 0.)
            } else if distance < 1. && point.y <= self.minimum + EPSILON {
                Vector::new(0., -1., 0.)
            } else {
                Vector::new(point.x, 0., point.z)
            }
        }

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
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
                (Point::new(1., 0., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0., 0., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0., 0., -5.), Vector::new(1., 1., 1.)),
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
                (Point::new(1., 0., -5.), Vector::new(0., 0., 1.), 5., 5.),
                (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 4., 6.),
                (
                    Point::new(0.5, 0., -5.),
                    Vector::new(0.1, 1., 1.),
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
                (Point::new(1., 0., 0.), Vector::new(1., 0., 0.)),
                (Point::new(0., 5., -1.), Vector::new(0., 0., -1.)),
                (Point::new(0., -2., 1.), Vector::new(0., 0., 1.)),
                (Point::new(-1., 1., 0.), Vector::new(-1., 0., 0.)),
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
                (Point::new(0., 1.5, 0.), Vector::new(0.1, 1., 0.), 0),
                (Point::new(0., 3., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 2., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 1., -5.), Vector::new(0., 0., 1.), 0),
                (Point::new(0., 1.5, -2.), Vector::new(0., 0., 1.), 2),
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
                (Point::new(0., 3., 0.), Vector::new(0., -1., 0.), 2),
                (Point::new(0., 3., -2.), Vector::new(0., -1., 2.), 2),
                (Point::new(0., 4., -2.), Vector::new(0., -1., 1.), 2),
                (Point::new(0., 0., -2.), Vector::new(0., 1., 2.), 2),
                (Point::new(0., -1., -2.), Vector::new(0., 1., 1.), 2),
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
                (Point::new(0., 1., 0.), Vector::new(0., -1., 0.)),
                (Point::new(0.5, 1., 0.), Vector::new(0., -1., 0.)),
                (Point::new(0., 1., 0.5), Vector::new(0., -1., 0.)),
                (Point::new(0., 2., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0.5, 2., 0.), Vector::new(0., 1., 0.)),
                (Point::new(0., 2., 0.5), Vector::new(0., 1., 0.)),
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
        parent: MaybeKeyRef,
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
            let parent = Arc::new(RwLock::new(parent));
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
                    vec![Intersection::new(t, Arc::new(self.clone()))]
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
                    intersections.push(Intersection::new(t0, Arc::new(self.clone())));
                }

                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection::new(t1, Arc::new(self.clone())));
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
                intersections.push(Intersection::new(t, Arc::new(self.clone())));
            }

            // check for an intersection with the upper end cap
            let t = (self.maximum.abs() - ray.origin.y) / ray.direction.y;
            if ray.check_cap(t, self.maximum.abs()) {
                intersections.push(Intersection::new(t, Arc::new(self.clone())));
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

        fn local_normal_at(&self, point: Point) -> Vector {
            let distance = point.x.powi(2) + point.z.powi(2);

            if distance < 1. && point.y >= self.maximum - EPSILON {
                Vector::new(0., 1., 0.)
            } else if distance < 1. && point.y <= self.minimum + EPSILON {
                Vector::new(0., -1., 0.)
            } else {
                let y = (point.x.powi(2) + point.z.powi(2)).sqrt();
                let y = if point.y > 0. { -y } else { y };
                Vector::new(point.x, y, point.z)
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

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
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
                (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 5., 5.),
                (
                    Point::new(0., 0., -5.),
                    Vector::new(1., 1., 1.),
                    8.66025,
                    8.66025,
                ),
                (
                    Point::new(1., 1., -5.),
                    Vector::new(-0.5, -1., 1.),
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
            let direction = Vector::new(0., 1., 1.).norm();
            let ray = Ray::new(Point::new(0., 0., -1.), direction);
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
                (Point::new(0., 0., -5.), Vector::new(0., 1., 0.), 0),
                (Point::new(0., 0., -0.25), Vector::new(0., 1., 1.), 2),
                (Point::new(0., 0., -0.25), Vector::new(0., 1., 0.), 4),
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
                (Point::new(0., 0., 0.), Vector::new(0., 0., 0.)),
                (
                    Point::new(1., 1., 1.),
                    Vector::new(1., -(2.0_f64.sqrt()), 1.),
                ),
                (Point::new(-1., -1., 0.), Vector::new(-1., 1., 0.)),
            ];

            for (point, expected_normal) in examples {
                let normal = shape.local_normal_at(point);
                assert_eq!(expected_normal, normal);
            }
        }
    }
}

pub mod triangle {
    use crate::EPSILON;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Triangle {
        id: usize,
        pub point1: Point,
        pub point2: Point,
        pub point3: Point,
        edge1: Vector,
        edge2: Vector,
        normal: Vector,
        material: Material,
        parent: MaybeKeyRef,
    }

    impl Triangle {
        pub fn new(point1: Point, point2: Point, point3: Point, material: Material) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);

            let edge1 = point2 - point1;
            let edge2 = point3 - point1;
            let normal = edge2.cross(edge1).norm();
            Self {
                id,
                point1,
                point2,
                point3,
                edge1,
                edge2,
                normal,
                material,
                parent: Arc::new(RwLock::new(None)),
            }
        }

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            let direction_cross_edge2 = ray.direction.cross(self.edge2);
            let determinant = self.edge1.dot(direction_cross_edge2);
            if determinant.abs() < EPSILON {
                return vec![];
            }

            let f = 1. / determinant;

            let point1_to_origin = ray.origin - self.point1;
            let u = f * point1_to_origin.dot(direction_cross_edge2);
            if !(0. ..=1.).contains(&u) {
                return vec![];
            }

            let origin_cross_edge1 = point1_to_origin.cross(self.edge1);
            let v = f * ray.direction.dot(origin_cross_edge1);
            if v < 0. || (u + v) > 1. {
                return vec![];
            }

            let t = f * self.edge2.dot(origin_cross_edge1);
            vec![Intersection::new(t, Arc::new(self.clone()))]
        }
    }

    impl Intersect for Triangle {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            Ok(self.local_intersect(ray))
        }
    }

    impl Shape for Triangle {
        fn id(&self) -> usize {
            self.id
        }

        fn local_normal_at(&self, _local_point: Point) -> Vector {
            self.normal
        }

        fn material(&self) -> Material {
            if let Some(parent) = self.parent() {
                // we know parent has to be a Group, so ok to unwrap here
                let group = parent.as_group().unwrap();
                if let Some(ref material) = group.material {
                    material.clone()
                } else {
                    self.material.clone()
                }
            } else {
                self.material.clone()
            }
        }

        fn set_material(&mut self, material: Material) {
            self.material = material;
        }

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn shape_type(&self) -> ShapeType {
            ShapeType::Triangle
        }

        fn transformation(&self) -> Matrix {
            Matrix::identity()
        }

        fn as_any(&self) -> &dyn Any {
            self
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn construct_triangle() {
            let p1 = Point::new(0., 1., 0.);
            let p2 = Point::new(-1., 0., 0.);
            let p3 = Point::new(1., 0., 0.);

            let triangle = Triangle::new(p1, p2, p3, Material::default());
            assert_eq!(triangle.point1, p1);
            assert_eq!(triangle.edge1, Vector::new(-1., -1., 0.));
            assert_eq!(triangle.edge2, Vector::new(1., -1., 0.));
            assert_eq!(triangle.normal, Vector::new(0., 0., -1.));
        }

        #[test]
        fn find_normal_on_triangle() {
            let triangle = Triangle::new(
                Point::new(0., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
                Material::default(),
            );
            let n1 = triangle.local_normal_at(Point::new(0., 0.5, 0.));
            let n2 = triangle.local_normal_at(Point::new(-0.5, 0.75, 0.));
            let n3 = triangle.local_normal_at(Point::new(0.5, 0.25, 0.));
            assert_eq!(n1, triangle.normal);
            assert_eq!(n2, triangle.normal);
            assert_eq!(n3, triangle.normal);
        }

        #[test]
        fn intersect_ray_parallel_to_triangle() {
            let triangle = Triangle::new(
                Point::new(0., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
                Material::default(),
            );
            let ray = Ray::new(Point::new(0., -1., -2.), Vector::new(0., 1., 0.));
            let xs = triangle.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_misses_p1_p3_edge() {
            let triangle = Triangle::new(
                Point::new(0., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
                Material::default(),
            );
            let ray = Ray::new(Point::new(1., 1., -2.), Vector::new(0., 0., 1.));
            let xs = triangle.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_misses_p1_p2_edge() {
            let triangle = Triangle::new(
                Point::new(0., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
                Material::default(),
            );

            let ray = Ray::new(Point::new(-1., 1., -2.), Vector::new(0., 0., 1.));
            let xs = triangle.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_misses_p2_p3_edge() {
            let triangle = Triangle::new(
                Point::new(0., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
                Material::default(),
            );

            let ray = Ray::new(Point::new(0., -1., -2.), Vector::new(0., 0., 1.));
            let xs = triangle.local_intersect(&ray);
            assert!(xs.is_empty());
        }

        #[test]
        fn ray_strikes_triangle() {
            let triangle = Triangle::new(
                Point::new(0., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
                Material::default(),
            );
            let ray = Ray::new(Point::new(0., 0.5, -2.), Vector::new(0., 0., 1.));
            let xs = triangle.local_intersect(&ray);
            assert_eq!(xs.len(), 1);
            assert_eq!(xs[0].t, 2.);
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
        parent: MaybeKeyRef,
        pub material: Option<Material>,
    }

    impl Group {
        pub fn new(
            transformation: Matrix,
            children: GroupChildren,
            parent: Option<DefaultKey>,
        ) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(1);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);

            let parent = Arc::new(RwLock::new(parent));
            let material = Some(Material::default());

            Self {
                id,
                transformation,
                children,
                parent,
                material,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }

        pub fn with_material(mut self, material: Material) -> Self {
            self.material = Some(material);
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

        /// Sets the parent of the child to the key of `self`, registers the child in `SHAPES`,
        /// and adds the key from `SHAPES` to `self`'s children.
        /// `self` must already be registered in `SHAPES` (either directly via `register_shape` or indirectly
        /// by being added as a child of another group).
        /// Panics if `self` is not in `SHAPES`.
        /// This is an O(n) operation because this uses `Shape.key()` underlyingly
        /// (see that method's documentation for an explanation of the time complexity).
        pub fn add_child(&self, child: &ShapeRef) {
            let self_key = self.key().expect("The group must be in shapes");
            child.set_parent(self_key);

            let child_key = register_shape(Arc::clone(child));
            let mut children = self.children.write();
            children.push(child_key);
        }

        pub fn get_child(&self, index: usize) -> Option<ShapeRef> {
            let children = self.children.read();
            let child_key = children.get(index)?;
            let shapes = SHAPES.read();
            let child = shapes.get(*child_key)?;
            Some(Arc::clone(child))
        }

        pub fn children(&self) -> Vec<ShapeRef> {
            let mut children_shapes = vec![];
            let children = self.children.read();
            for child_index in 0..children.len() {
                let child = self.get_child(child_index).unwrap();
                children_shapes.push(child);
            }
            children_shapes
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

        fn local_normal_at(&self, _: Point) -> Vector {
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

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
        }
    }

    impl Intersect for Group {
        fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
            let ray = ray.transform(self.transformation.inverse()?);
            self.local_intersect(&ray)
        }
    }

    impl PartialEq for Group {
        fn eq(&self, other: &Self) -> bool {
            self.id == other.id
        }
    }

    impl Eq for Group {}

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
            register_shape(group_clone);

            let shape: ShapeRef = Arc::new(shape);
            group.add_child(&shape);

            let children = group.children.read();
            assert!(!children.is_empty());

            let child = group.get_child(0).unwrap();

            assert_eq!(shape.id(), child.id());
            assert_eq!(shape.shape_type(), child.shape_type());

            let shapes = SHAPES.read();
            let parent_key = child.parent_key().unwrap();
            let parent = shapes.get(parent_key).unwrap();
            assert_eq!(parent.id(), group.id);
        }

        #[test]
        fn intersect_ray_with_empty_group() {
            let group = Group::default();
            let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
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
            register_shape(group_clone);

            let s1: ShapeRef = Arc::new(s1);
            let s2: ShapeRef = Arc::new(s2);
            let s3: ShapeRef = Arc::new(s3);

            group.add_child(&s1);
            group.add_child(&s2);
            group.add_child(&s3);

            let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
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
            register_shape(group_clone);

            let shape: ShapeRef = Arc::new(shape);
            group.add_child(&shape);

            let ray = Ray::new(Point::new(10., 0., -10.), Vector::new(0., 0., 1.));
            let xs = group.intersect(&ray).unwrap();
            assert_eq!(xs.len(), 2);
        }
    }
}

#[cfg(test)]
mod tests {

    use std::{
        f64::consts::{FRAC_1_SQRT_2, PI},
        thread::sleep,
        time::Duration,
    };

    use super::{group::Group, sphere::Sphere, *};

    #[derive(Debug, Clone)]
    pub struct TestShape {
        transformation: Matrix,
        material: Material,
        parent: MaybeKeyRef,
    }

    impl TestShape {
        pub fn new(transformation: Matrix, material: Material, parent: Option<DefaultKey>) -> Self {
            let parent = Arc::new(RwLock::new(parent));
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

        fn local_normal_at(&self, local_point: Point) -> Vector {
            Vector::new(local_point.x, local_point.y, local_point.z)
        }

        fn parent_key(&self) -> Option<DefaultKey> {
            *self.parent.read()
        }

        fn set_parent(&self, parent: DefaultKey) {
            *self.parent.write() = Some(parent);
        }

        fn as_any(&self) -> &dyn Any {
            self
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
            .normal_at(Point::new(0., 1.70711, -FRAC_1_SQRT_2))
            .unwrap();
        assert_eq!(normal, Vector::new(0., FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn compute_normal_on_transformed_shape() {
        let transformation = Matrix::identity().rotate_z(PI / 2.).scale(1., 0.5, 1.);
        let shape = TestShape::default().with_transformation(transformation);
        let val = 2.0_f64.sqrt() / 2.;
        let normal = shape.normal_at(Point::new(0., val, -val)).unwrap();
        assert_eq!(normal, Vector::new(0., 0.97014, -0.24254));
    }

    #[test]
    fn shape_has_parent_field() {
        let shape = TestShape::default();
        assert!(shape.parent_key().is_none());
    }

    #[test]
    fn convert_point_from_world_to_object_space() {
        while SHAPES.is_locked() {
            sleep(Duration::from_millis(100));
        }

        let sphere: ShapeRef =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(5., 0., 0.)));

        let group1 = Arc::new(Group::default().with_transformation(Matrix::rotation_y(PI / 2.)));
        let group2: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::scaling(2., 2., 2.)));

        let group1_clone: ShapeRef = Arc::clone(&group1) as ShapeRef;
        register_shape(group1_clone);

        group1.add_child(&group2);

        let group2 = group2.as_group().unwrap();
        group2.add_child(&sphere);

        let point = sphere.world_to_object(Point::new(-2., 0., -10.)).unwrap();
        assert_eq!(point, Point::new(0., 0., -1.));
    }

    #[test]
    fn convert_normal_from_object_to_world_space() {
        while SHAPES.is_locked() {
            sleep(Duration::from_millis(100));
        }

        let group1 = Arc::new(Group::default().with_transformation(Matrix::rotation_y(PI / 2.)));
        let group1_clone = Arc::clone(&group1) as ShapeRef;
        register_shape(group1_clone);

        let group2: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::scaling(1., 2., 3.)));
        group1.add_child(&group2);

        let sphere: ShapeRef =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(5., 0., 0.)));
        let group2 = group2.as_group().unwrap();
        group2.add_child(&sphere);

        let val = 3.0_f64.sqrt() / 3.;
        let normal = sphere.normal_to_world(Vector::new(val, val, val)).unwrap();

        // need more precise values (5 significant figures)
        // than the book provides
        // cf. book values (0.2857, 0.4286, -0.8571)
        assert_eq!(normal, Vector::new(0.28571, 0.42857, -0.85714));
    }

    #[test]
    fn find_normal_on_child_object() {
        while SHAPES.is_locked() {
            sleep(Duration::from_millis(100));
        }

        let group1 = Arc::new(Group::default().with_transformation(Matrix::rotation_y(PI / 2.)));
        register_shape(Arc::clone(&group1) as ShapeRef);

        let group2: ShapeRef =
            Arc::new(Group::default().with_transformation(Matrix::scaling(1., 2., 3.)));
        group1.add_child(&group2);

        let sphere: ShapeRef =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(5., 0., 0.)));
        let group2 = group2.as_group().unwrap();
        group2.add_child(&sphere);

        let normal = sphere
            .normal_at(Point::new(1.7321, 1.1547, -5.5774))
            .unwrap();
        // cf. book values (0.2857, 0.4286, -0.8571)
        assert_eq!(normal, Vector::new(0.2857, 0.42854, -0.85716));
    }

    #[test]
    fn get_back_original_type() {
        let shape: Arc<dyn Shape> = Arc::new(TestShape::default());
        let _concrete_shape = shape.as_any().downcast_ref::<TestShape>().unwrap();
    }
}
