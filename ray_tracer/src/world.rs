use std::sync::Arc;

use crate::{
    intersection::{hit, Computations, Intersect, Intersection},
    light::{Light, PointLight},
    material::Material,
    shape::sphere::Sphere,
    shape::ShapeRef,
    Color, Matrix, Point, Ray, Result,
};

#[derive(Debug)]
pub struct World {
    pub objects: Vec<ShapeRef>,
    pub lights: Vec<Box<dyn Light>>,
}

impl World {
    pub fn new(objects: Vec<ShapeRef>, lights: Vec<Box<dyn Light>>) -> Self {
        Self { objects, lights }
    }

    /// Constructs an empty world with no objects and no light.
    pub fn new_empty() -> Self {
        Self::new(vec![], vec![])
    }

    /// Intersects the world with the given ray and returns
    /// the color at the resulting intersection.
    pub fn color_at(&self, ray: &Ray, remaining: usize) -> Result<Color> {
        let xs = self.intersect(ray)?;
        if let Some(hit) = hit(&xs) {
            let comps = hit.prepare_computations(ray, &xs)?;
            Ok(self.shade_hit(&comps, remaining)?)
        } else {
            Ok(Color::default())
        }
    }

    /// Returns the color at the intersection encapsulated by comps,
    /// in the given world. Returns black if the world has no light source.
    fn shade_hit(&self, comps: &Computations, remaining: usize) -> Result<Color> {
        let mut colors = vec![];

        for light in &self.lights {
            let light_intensity = light.intensity_at(comps.over_point, self)?;
            let surface_color = comps.object.material().lighting(
                comps.object,
                &**light,
                comps.over_point,
                comps.eye_vector,
                comps.normal_vector,
                light_intensity,
            )?;

            let reflected_color = self.reflected_color(comps, remaining)?;
            let refracted_color = self.refracted_color(comps, remaining)?;

            let material = comps.object.material();
            let color = if material.reflectivity > 0. && material.transparency > 0. {
                let reflectance = comps.schlick();
                surface_color + reflected_color * reflectance + refracted_color * (1. - reflectance)
            } else {
                surface_color + reflected_color + refracted_color
            };

            colors.push(color);
        }

        if let Some(final_color) = colors.into_iter().reduce(|acc, element| acc + element) {
            Ok(final_color)
        } else {
            // there are no lights, so everything is black
            Ok(Color::black())
        }
    }

    pub fn is_shadowed(&self, point: Point, light_position: Point) -> Result<bool> {
        let distance_vector = light_position - point;
        let distance = distance_vector.magnitude();
        let direction = distance_vector.norm();

        let ray = Ray::new(point, direction);
        let intersections = self.intersect(&ray)?;

        if let Some(hit) = hit(&intersections) {
            Ok(hit.t < distance && hit.object.material().casts_shadow)
        } else {
            Ok(false)
        }
    }

    fn reflected_color(&self, comps: &Computations, remaining: usize) -> Result<Color> {
        if remaining == 0 {
            return Ok(Color::black());
        }

        let material_reflective = comps.object.material().reflectivity;
        if material_reflective == 0. {
            return Ok(Color::black());
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflect_vector);
        let color = self.color_at(&reflect_ray, remaining - 1)?;
        Ok(color * material_reflective)
    }

    fn refracted_color(&self, comps: &Computations, remaining: usize) -> Result<Color> {
        if comps.object.material().transparency == 0. || remaining == 0 {
            return Ok(Color::black());
        }

        // find the ratio of first index of refraction to the second
        // (this is inverted from the definition of Snell's law)
        if let (Some(n1), Some(n2)) = (comps.n1, comps.n2) {
            let n_ratio = n1 / n2;

            // cos(theta_i) is the same as the dot product of the two vectors
            let cos_i = comps.eye_vector.dot(comps.normal_vector);

            // find sin(theta_t)^2 via trigonometric identity
            let sin2_t = n_ratio.powi(2) * (1. - cos_i.powi(2));

            if sin2_t > 1. {
                return Ok(Color::black());
            }

            // find cos(theta_t) via trigonometric identity
            let cos_t = (1. - sin2_t).sqrt();

            // compute the direction of the refracted ray
            let direction =
                comps.normal_vector * (n_ratio * cos_i - cos_t) - comps.eye_vector * n_ratio;

            // create the refracted ray
            let refracted_ray = Ray::new(comps.under_point, direction);

            // find the color of the refracted ray,
            // making sure to multiply by the transparency value
            // to account fo any opacity
            let color =
                self.color_at(&refracted_ray, remaining)? * comps.object.material().transparency;
            return Ok(color);
        }

        Ok(Color::white())
    }
}

impl Default for World {
    /// Constructs the default world with a light source at (-10, 10, -10)
    /// and two concentric spheres, where the outermost is a unit sphere
    /// and the inntermost has a radius of 0.5. Both lie at the origin.
    fn default() -> Self {
        let light = Box::new(PointLight::new(Point::new(-10., 10., -10.), Color::white()));

        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Arc::new(Sphere::default().with_material(material));
        let sphere2 = Arc::new(Sphere::default().with_transformation(transformation));

        Self {
            lights: vec![light],
            objects: vec![sphere1, sphere2],
        }
    }
}

impl Intersect for World {
    fn intersect(&self, ray: &Ray) -> Result<Vec<Intersection>> {
        let mut xs = vec![];
        for object in &self.objects {
            let mut intersections = object.intersect(ray)?;
            xs.append(&mut intersections);
        }
        xs.sort_unstable_by(|a, b| a.t.total_cmp(&b.t));
        Ok(xs)
    }
}

#[cfg(test)]
mod tests {
    use crate::{pattern, shape::plane::Plane, Vector};

    use super::*;

    #[test]
    fn create_a_world() {
        let world = World::new_empty();
        assert!(world.objects.is_empty());
        assert_eq!(world.lights, vec![]);
    }

    #[test]
    fn create_default_world() {
        let mut world = World::default();
        let light = Box::new(PointLight::new(Point::new(-10., 10., -10.), Color::white()))
            as Box<dyn Light>;
        assert_eq!(world.lights, vec![light]);

        let color = Color::new(0.8, 1.0, 0.6);
        let diffuse = 0.7;
        let specular = 0.2;
        let material = Material {
            color,
            diffuse,
            specular,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);
        let sphere1 = Sphere::default().with_material(material);
        let sphere2 = Sphere::default().with_transformation(transformation);

        // kind of an ugly hack to deal with the fact
        // that we can't easily clone trait objects
        // and that our Sphere values would get moved as soon as we pass them to the world
        // also, we need to explicitly assign those objects to the sphere
        // because if we just rely on the default ones, they'd have different id's
        let boxed_sphere1: ShapeRef = Arc::new(sphere1.clone());
        let boxed_sphere2: ShapeRef = Arc::new(sphere2.clone());
        world.objects = vec![boxed_sphere1, boxed_sphere2];

        let boxed_sphere1: ShapeRef = Arc::new(sphere1);
        let boxed_sphere2: ShapeRef = Arc::new(sphere2);
        assert!(world.objects.contains(&boxed_sphere1));
        assert!(world.objects.contains(&boxed_sphere2));
    }

    #[test]
    fn intersect_world_with_ray() {
        let world = World::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = world.intersect(&ray).unwrap();
        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0].t, 4.);
        assert_eq!(xs[1].t, 4.5);
        assert_eq!(xs[2].t, 5.5);
        assert_eq!(xs[3].t, 6.);
    }

    #[test]
    fn shading_intersection() {
        let world = World::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = &world.objects[0];
        let intersection = Intersection::new(4., shape.clone());
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        let color = world.shade_hit(&comps, 5).unwrap();
        assert_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_intersection_from_inside() {
        let world = World {
            lights: vec![Box::new(PointLight::new(
                Point::new(0., 0.25, 0.),
                Color::white(),
            ))],
            ..Default::default()
        };
        let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));

        let shape = &world.objects[1];
        let intersection = Intersection::new(0.5, shape.clone());
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        let color = world.shade_hit(&comps, 5).unwrap();

        let val = 0.90498;
        assert_eq!(color, Color::new(val, val, val));
    }

    #[test]
    fn shade_hit_is_given_intersection_in_shadow() {
        let light = Box::new(PointLight::new(Point::new(0., 0., -10.), Color::white()));
        let sphere1 = Arc::<Sphere>::default();
        let sphere2 =
            Arc::new(Sphere::default().with_transformation(Matrix::translation(0., 0., 10.)));
        let world = World {
            lights: vec![light],
            objects: vec![sphere1, sphere2.clone()],
        };

        let ray = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let intersection = Intersection::new(4., sphere2);
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        let color = world.shade_hit(&comps, 5).unwrap();
        assert_eq!(color, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn color_when_ray_misses() {
        // when the ray fails to intersect anything, the returned color should be black
        let world = World::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        let color = world.color_at(&ray, 5).unwrap();
        assert_eq!(color, Color::default());
    }

    #[test]
    fn color_when_ray_hits() {
        // our ray intersects the outermost sphere in the default world
        let world = World::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let color = world.color_at(&ray, 5).unwrap();
        assert_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn color_with_intersection_behind_ray() {
        // color_at should use `hit` when computing the color
        // We put the ray inside the outer sphere, but outside the inner sphere,
        // pointing at the inner sphere.
        // We expect the hit to be on the inner sphere, and so return its color
        let mut world = World::default();
        let outer = Arc::get_mut(&mut world.objects[0]).unwrap();
        let mut outer_material = outer.material();
        outer_material.ambient = 1.;
        outer.set_material(outer_material);

        let inner = Arc::get_mut(&mut world.objects[1]).unwrap();
        let mut inner_material = inner.material();
        inner_material.ambient = 1.;
        inner.set_material(inner_material.clone());

        let ray = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let color = world.color_at(&ray, 5).unwrap();
        assert_eq!(color, inner_material.color);
    }

    #[test]
    fn reflected_color_for_nonreflective_material() {
        // we have to replicate all this code here to avoid a borrow checker error
        // that would result from a simultaneous mutable and immutable borrow of world
        // which would happen if we went the route of constructing a default world
        // then getting a mutable reference to its second shape
        // and trying to set the material on that directly
        let light = Box::new(PointLight::new(Point::new(-10., 10., -10.), Color::white()));

        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ambient: 1.,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Arc::new(Sphere::default().with_material(material));
        let sphere2 = Arc::new(Sphere::default().with_transformation(transformation));

        let world = World {
            lights: vec![light],
            objects: vec![sphere1, sphere2],
        };
        let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));

        let intersection = Intersection::new(1., world.objects[1].clone());
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();
        let color = world.reflected_color(&comps, 5).unwrap();
        assert_eq!(color, Color::black());
    }

    #[test]
    fn reflected_color_for_reflective_material() {
        // ditto about repeating the code to prevent running afoul of the borrow checker
        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Arc::new(Sphere::default().with_material(material));
        let sphere2 = Arc::new(Sphere::default().with_transformation(transformation));

        let shape = Plane::default()
            .with_material(Material {
                reflectivity: 0.5,
                ..Default::default()
            })
            .with_transformation(Matrix::translation(0., -1., 0.));
        let world = World {
            objects: vec![sphere1, sphere2, Arc::new(shape)],
            ..Default::default()
        };

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Point::new(0., 0., -3.), Vector::new(0., -val, val));
        let intersection = Intersection::new(2.0_f64.sqrt(), world.objects[2].clone());
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();

        let color = world.reflected_color(&comps, 5).unwrap();

        // the expected values are slightly modified from the book
        // to account for floating point errors
        // which occur because the values given by the book are off from the ones we get here
        // by more than EPSILON
        // cf. book values: color(0.19032, 0.2379, 0.14274)
        assert_eq!(color, Color::new(0.19033, 0.23792, 0.14274));
    }

    #[test]
    fn shade_hit_with_reflective_material() {
        // ditto about repeating the code to prevent running afoul of the borrow checker
        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Arc::new(Sphere::default().with_material(material));
        let sphere2 = Arc::new(Sphere::default().with_transformation(transformation));

        let shape = Plane::default()
            .with_material(Material {
                reflectivity: 0.5,
                ..Default::default()
            })
            .with_transformation(Matrix::translation(0., -1., 0.));
        let world = World {
            objects: vec![sphere1, sphere2, Arc::new(shape)],
            ..Default::default()
        };

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Point::new(0., 0., -3.), Vector::new(0., -val, val));
        let intersection = Intersection::new(2.0_f64.sqrt(), world.objects[2].clone());
        let comps = intersection
            .prepare_computations(&ray, &[intersection.to_owned()])
            .unwrap();

        let color = world.shade_hit(&comps, 5).unwrap();

        // ditto as in the above test about using modified values
        // cf. book values: color(0.87677, 0.92436, 0.82918)
        assert_eq!(color, Color::new(0.87675, 0.92434, 0.82917));
    }

    #[test]
    fn color_at_with_mutually_reflective_surfaces() {
        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Arc::new(Sphere::default().with_material(material));
        let sphere2 = Arc::new(Sphere::default().with_transformation(transformation));

        let light = Box::new(PointLight::new(
            Point::new(0., 0., 0.),
            Color::new(1., 1., 1.),
        ));

        let lower = Arc::new(
            Plane::default()
                .with_material(Material {
                    reflectivity: 1.,
                    ..Default::default()
                })
                .with_transformation(Matrix::translation(0., -1., 0.)),
        );
        let upper = Arc::new(
            Plane::default()
                .with_material(Material {
                    reflectivity: 1.,
                    ..Default::default()
                })
                .with_transformation(Matrix::translation(0., 1., 0.)),
        );

        let world = World::new(vec![sphere1, sphere2, lower, upper], vec![light]);

        let ray = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));

        // should terminate successfully without entering infinite recursion
        let _ = world.color_at(&ray, 5).unwrap();
    }

    #[test]
    fn reflected_color_at_maximum_recursive_depth() {
        let mut world = World::default();
        let shape = Arc::new(
            Plane::default()
                .with_material(Material {
                    reflectivity: 0.5,
                    ..Default::default()
                })
                .with_transformation(Matrix::translation(0., -1., 0.)),
        );
        world.objects.push(shape.clone());

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Point::new(0., 0., -3.), Vector::new(0., -val, val));
        let intersection = Intersection::new(2.0_f64.sqrt(), shape);

        let comps = intersection
            .prepare_computations(&ray, &[intersection.clone()])
            .unwrap();
        let color = world.reflected_color(&comps, 0).unwrap();
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_with_opaque_surface() {
        let world = World::default();
        let shape = &world.objects[0];
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = &[
            Intersection::new(4., shape.clone()),
            Intersection::new(6., shape.clone()),
        ];
        let comps = xs[0].prepare_computations(&ray, xs).unwrap();
        let color = world.refracted_color(&comps, 5).unwrap();
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_at_maximum_recursive_depth() {
        let mut world = World::default();
        let shape = Arc::get_mut(&mut world.objects[0]).unwrap();
        let mut material = shape.material();
        material.transparency = 1.;
        material.refractive_index = 1.5;
        shape.set_material(material);

        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Arc::clone(&world.objects[0]);
        let xs = &[
            Intersection::new(4., shape.clone()),
            Intersection::new(6., shape.clone()),
        ];

        let comps = xs[0].prepare_computations(&ray, xs).unwrap();
        let color = world.refracted_color(&comps, 0).unwrap();
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_under_total_internal_reflection() {
        let mut world = World::default();
        let shape = Arc::get_mut(&mut world.objects[0]).unwrap();
        let mut material = shape.material();
        material.transparency = 1.;
        material.refractive_index = 1.5;
        shape.set_material(material);

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Point::new(0., 0., val), Vector::new(0., 1., 0.));
        let shape = Arc::clone(&world.objects[0]);
        let xs = &[
            Intersection::new(-val, shape.clone()),
            Intersection::new(val, shape.clone()),
        ];

        // this time we're inside the sphere,
        // so we need to look at the second intersection
        let comps = xs[1].prepare_computations(&ray, xs).unwrap();
        let color = world.refracted_color(&comps, 5).unwrap();
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_with_refracted_ray() {
        let mut world = World::default();

        // we have to do this in separate blocks
        // to work around the ban on multiple mutable references to world.objects
        {
            let a = Arc::get_mut(&mut world.objects[0]).unwrap();
            let mut a_material = a.material();
            a_material.ambient = 1.;
            a_material.pattern = Some(Box::<pattern::tests::TestPattern>::default());
            a.set_material(a_material);
        }

        {
            let b = Arc::get_mut(&mut world.objects[1]).unwrap();
            let mut b_material = b.material();
            b_material.transparency = 1.;
            b_material.refractive_index = 1.5;
            b.set_material(b_material);
        }

        let a = &world.objects[0];
        let b = &world.objects[1];

        let ray = Ray::new(Point::new(0., 0., 0.1), Vector::new(0., 1., 0.));
        let xs = &[
            Intersection::new(-0.9899, a.clone()),
            Intersection::new(-0.4899, b.clone()),
            Intersection::new(0.4899, b.clone()),
            Intersection::new(0.9899, a.clone()),
        ];

        let comps = xs[2].prepare_computations(&ray, xs).unwrap();
        let color = world.refracted_color(&comps, 5).unwrap();

        // cf. book values color(0, 0.99888, 0.04725)
        assert_eq!(color, Color::new(0., 0.99887, 0.04722));
    }

    #[test]
    fn shade_hit_with_transparent_material() {
        let mut world = World::default();

        let floor = Plane::default()
            .with_transformation(Matrix::translation(0., -1., 0.))
            .with_material(Material {
                transparency: 0.5,
                refractive_index: 1.5,
                ..Default::default()
            });
        world.objects.push(Arc::new(floor.clone()));

        let ball = Sphere::new(
            Matrix::translation(0., -3.5, -0.5),
            Material {
                color: Color::new(1., 0., 0.),
                ambient: 0.5,
                ..Default::default()
            },
            None,
        );
        world.objects.push(Arc::new(ball));

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Point::new(0., 0., -3.), Vector::new(0., -val, val));
        let xs = &[Intersection::new(2.0_f64.sqrt(), Arc::new(floor))];
        let comps = xs[0].prepare_computations(&ray, xs).unwrap();
        let color = world.shade_hit(&comps, 5).unwrap();
        assert_eq!(color, Color::new(0.93642, 0.68642, 0.68642));
    }

    #[test]
    fn shade_hit_with_reflective_transparent_material() {
        let mut world = World::default();

        let floor = Plane::default()
            .with_transformation(Matrix::translation(0., -1., 0.))
            .with_material(Material {
                reflectivity: 0.5,
                transparency: 0.5,
                refractive_index: 1.5,
                ..Default::default()
            });
        world.objects.push(Arc::new(floor.clone()));

        let ball = Sphere::new(
            Matrix::translation(0., -3.5, -0.5),
            Material {
                color: Color::new(1., 0., 0.),
                ambient: 0.5,
                ..Default::default()
            },
            None,
        );
        world.objects.push(Arc::new(ball));

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Point::new(0., 0., -3.), Vector::new(0., -val, val));

        let xs = &[Intersection::new(2.0_f64.sqrt(), Arc::new(floor))];
        let comps = xs[0].prepare_computations(&ray, xs).unwrap();
        let color = world.shade_hit(&comps, 5).unwrap();

        assert_eq!(color, Color::new(0.93391, 0.69643, 0.69243));
    }

    #[test]
    fn is_shadow_tests_for_occlusion_between_2_points() {
        let world = World::default();
        let light_position = Point::new(-10., -10., -10.);
        let examples = [
            (Point::new(-10., -10., 10.), false),
            (Point::new(10., 10., 10.), true),
            (Point::new(-20., -20., -20.), false),
            (Point::new(-5., -5., -5.), false),
        ];

        for (point, expected) in examples {
            let result = world.is_shadowed(point, light_position).unwrap();
            assert_eq!(result, expected);
        }
    }
}
