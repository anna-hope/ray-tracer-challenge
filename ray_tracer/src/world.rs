use crate::{
    intersection::{hit, Computations, Intersect, Intersection},
    light::PointLight,
    material::Material,
    shape::sphere::Sphere,
    shape::Shape,
    Color, Matrix, Ray, Result, Tuple,
};

#[derive(Debug)]
pub struct World {
    pub objects: Vec<Box<dyn Shape>>,
    pub light: Option<PointLight>,
}

impl World {
    /// Constructs an empty world with no objects and no light.
    pub fn new(objects: Vec<Box<dyn Shape>>, light: Option<PointLight>) -> Self {
        Self { objects, light }
    }

    pub fn new_empty() -> Self {
        Self::new(vec![], None)
    }

    /// Returns the color at the intersection encapsulated by comps,
    /// in the given world. Returns black if the world has no light source.
    fn shade_hit(&self, comps: &Computations) -> Result<Color> {
        // possible future improvement: support multiple light sources
        // (p. 96 of the book)
        if let Some(light) = self.light {
            let in_shadow = self.is_shadowed(comps.over_point)?;
            let surface_color = comps.object.material().lighting(
                comps.object,
                light,
                comps.over_point,
                comps.eye_vector,
                comps.normal_vector,
                in_shadow,
            )?;

            let reflected_color = self.reflected_color(comps)?;

            Ok(surface_color + reflected_color)
        } else {
            Ok(Color::default())
        }
    }

    /// Intersects the world with the given ray and returns
    /// the color at the resulting intersection.
    pub fn color_at(&self, ray: &Ray) -> Result<Color> {
        let xs = self.intersect(ray)?;
        if let Some(hit) = hit(&xs) {
            let comps = hit.prepare_computations(ray)?;
            Ok(self.shade_hit(&comps)?)
        } else {
            Ok(Color::default())
        }
    }

    fn is_shadowed(&self, point: Tuple) -> Result<bool> {
        if let Some(light) = self.light {
            let distance_vector = light.position - point;
            let distance = distance_vector.magnitude();
            let direction = distance_vector.norm();

            let ray = Ray::new(point, direction);
            let intersections = self.intersect(&ray)?;

            if let Some(hit) = hit(&intersections) {
                Ok(hit.t < distance)
            } else {
                Ok(false)
            }
        } else {
            // if there is no light, everything is shadowed
            Ok(true)
        }
    }

    fn reflected_color(&self, comps: &Computations) -> Result<Color> {
        let material_reflective = comps.object.material().reflective;
        if material_reflective == 0. {
            return Ok(Color::black());
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflect_vector);
        let color = self.color_at(&reflect_ray)?;
        Ok(color * material_reflective)
    }
}

impl Default for World {
    /// Constructs the default world with a light source at (-10, 10, -10)
    /// and two concentric spheres, where the outermost is a unit sphere
    /// and the inntermost has a radius of 0.5. Both lie at the origin.
    fn default() -> Self {
        let light = PointLight::new(Tuple::point(-10., 10., -10.), Color::white());

        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Box::new(Sphere::new().with_material(material));
        let sphere2 = Box::new(Sphere::new().with_transformation(transformation));

        Self {
            light: Some(light),
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
    use crate::shape::plane::Plane;

    use super::*;

    #[test]
    fn create_a_world() {
        let world = World::new_empty();
        assert!(world.objects.is_empty());
        assert_eq!(world.light, None);
    }

    #[test]
    fn create_default_world() {
        let mut world = World::default();
        let light = PointLight::new(Tuple::point(-10., 10., -10.), Color::white());
        assert_eq!(world.light, Some(light));

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
        let sphere1 = Sphere::new().with_material(material);
        let sphere2 = Sphere::new().with_transformation(transformation);

        // kind of an ugly hack to deal with the fact
        // that we can't easily clone trait objects
        // and that our Sphere values would get moved as soon as we pass them to the world
        // also, we need to explicitly assign those objects to the sphere
        // because if we just rely on the default ones, they'd have different id's
        let boxed_sphere1: Box<dyn Shape> = Box::new(sphere1.clone());
        let boxed_sphere2: Box<dyn Shape> = Box::new(sphere2.clone());
        world.objects = vec![boxed_sphere1, boxed_sphere2];

        let boxed_sphere1: Box<dyn Shape> = Box::new(sphere1.clone());
        let boxed_sphere2: Box<dyn Shape> = Box::new(sphere2.clone());
        assert!(world.objects.contains(&boxed_sphere1));
        assert!(world.objects.contains(&boxed_sphere2));
    }

    #[test]
    fn intersect_world_with_ray() {
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
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
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let shape = &world.objects[0];
        let intersection = shape.arbitrary_intersection(4.);
        let comps = intersection.prepare_computations(&ray).unwrap();
        let color = world.shade_hit(&comps).unwrap();
        assert_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_intersection_from_inside() {
        let world = World {
            light: Some(PointLight::new(Tuple::point(0., 0.25, 0.), Color::white())),
            ..Default::default()
        };
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));

        let shape = &world.objects[1];
        let intersection = shape.arbitrary_intersection(0.5);
        let comps = intersection.prepare_computations(&ray).unwrap();
        let color = world.shade_hit(&comps).unwrap();

        let val = 0.90498;
        assert_eq!(color, Color::new(val, val, val));
    }

    #[test]
    fn shade_hit_is_given_intersection_in_shadow() {
        let light = PointLight::new(Tuple::point(0., 0., -10.), Color::white());
        let sphere1 = Sphere::new();
        let sphere2 = Sphere::new().with_transformation(Matrix::translation(0., 0., 10.));
        let world = World {
            light: Some(light),
            objects: vec![Box::new(sphere1), Box::new(sphere2.clone())],
        };

        let ray = Ray::new(Tuple::point(0., 0., 5.), Tuple::vector(0., 0., 1.));
        let intersection = Intersection::new(4., &sphere2);
        let comps = intersection.prepare_computations(&ray).unwrap();
        let color = world.shade_hit(&comps).unwrap();
        assert_eq!(color, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn color_when_ray_misses() {
        // when the ray fails to intersect anything, the returned color should be black
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 1., 0.));
        let color = world.color_at(&ray).unwrap();
        assert_eq!(color, Color::default());
    }

    #[test]
    fn color_when_ray_hits() {
        // our ray intersects the outermost sphere in the default world
        let world = World::default();
        let ray = Ray::new(Tuple::point(0., 0., -5.), Tuple::vector(0., 0., 1.));
        let color = world.color_at(&ray).unwrap();
        assert_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn color_with_intersection_behind_ray() {
        // color_at should use `hit` when computing the color
        // We put the ray inside the outer sphere, but outside the inner sphere,
        // pointing at the inner sphere.
        // We expect the hit to be on the inner sphere, and so return its color
        let mut world = World::default();
        let outer = &mut world.objects[0];
        let mut outer_material = outer.material();
        outer_material.ambient = 1.;
        outer.set_material(outer_material);

        let inner = &mut world.objects[1];
        let mut inner_material = inner.material();
        inner_material.ambient = 1.;
        inner.set_material(inner_material.clone());

        let ray = Ray::new(Tuple::point(0., 0., 0.75), Tuple::vector(0., 0., -1.));
        let color = world.color_at(&ray).unwrap();
        assert_eq!(color, inner_material.color);
    }

    #[test]
    fn no_shadow_when_nothing_is_collinear_with_point_and_light() {
        let world = World::default();
        let point = Tuple::point(0., 10., 0.);
        let is_shadowed = world.is_shadowed(point).unwrap();
        assert!(!is_shadowed);
    }

    #[test]
    fn shadow_when_object_is_between_point_and_light() {
        let world = World::default();
        let point = Tuple::point(10., -10., 10.);
        assert!(world.is_shadowed(point).unwrap());
    }

    #[test]
    fn no_shadow_when_object_is_behind_light() {
        let world = World::default();
        let point = Tuple::point(-20., 20., -20.);
        assert!(!world.is_shadowed(point).unwrap());
    }

    #[test]
    fn no_shadow_when_object_is_behind_point() {
        let world = World::default();
        let point = Tuple::point(-2., -2., -2.);
        assert!(!world.is_shadowed(point).unwrap());
    }

    #[test]
    fn reflected_color_for_nonreflective_material() {
        // we have to replicate all this code here to avoid a borrow checker error
        // that would result from a simultaneous mutable and immutable borrow of world
        // which would happen if we went the route of constructing a default world
        // then getting a mutable reference to its second shape
        // and trying to set the material on that directly
        let light = PointLight::new(Tuple::point(-10., 10., -10.), Color::white());

        let material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ambient: 1.,
            ..Default::default()
        };

        let transformation = Matrix::scaling(0.5, 0.5, 0.5);

        let sphere1 = Box::new(Sphere::new().with_material(material));
        let sphere2 = Box::new(Sphere::new().with_transformation(transformation));

        let world = World {
            light: Some(light),
            objects: vec![sphere1, sphere2],
        };
        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 0., 1.));

        let intersection = &world.objects[1].arbitrary_intersection(1.);
        let comps = intersection.prepare_computations(&ray).unwrap();
        let color = world.reflected_color(&comps).unwrap();
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

        let sphere1 = Box::new(Sphere::new().with_material(material));
        let sphere2 = Box::new(Sphere::new().with_transformation(transformation));

        let shape = Plane::new()
            .with_material(Material {
                reflective: 0.5,
                ..Default::default()
            })
            .with_transformation(Matrix::translation(0., -1., 0.));
        let world = World {
            objects: vec![sphere1, sphere2, Box::new(shape)],
            ..Default::default()
        };

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Tuple::point(0., 0., -3.), Tuple::vector(0., -val, val));
        let intersection = &world.objects[2].arbitrary_intersection(2.0_f64.sqrt());
        let comps = intersection.prepare_computations(&ray).unwrap();

        let color = world.reflected_color(&comps).unwrap();

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

        let sphere1 = Box::new(Sphere::new().with_material(material));
        let sphere2 = Box::new(Sphere::new().with_transformation(transformation));

        let shape = Plane::new()
            .with_material(Material {
                reflective: 0.5,
                ..Default::default()
            })
            .with_transformation(Matrix::translation(0., -1., 0.));
        let world = World {
            objects: vec![sphere1, sphere2, Box::new(shape)],
            ..Default::default()
        };

        let val = 2.0_f64.sqrt() / 2.;
        let ray = Ray::new(Tuple::point(0., 0., -3.), Tuple::vector(0., -val, val));
        let intersection = &world.objects[2].arbitrary_intersection(2.0_f64.sqrt());
        let comps = intersection.prepare_computations(&ray).unwrap();

        let color = world.shade_hit(&comps).unwrap();

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

        let sphere1 = Box::new(Sphere::new().with_material(material));
        let sphere2 = Box::new(Sphere::new().with_transformation(transformation));

        let light = PointLight::new(Tuple::point(0., 0., 0.), Color::new(1., 1., 1.));

        let lower = Box::new(
            Plane::new()
                .with_material(Material {
                    reflective: 1.,
                    ..Default::default()
                })
                .with_transformation(Matrix::translation(0., -1., 0.)),
        );
        let upper = Box::new(
            Plane::new()
                .with_material(Material {
                    reflective: 1.,
                    ..Default::default()
                })
                .with_transformation(Matrix::translation(0., 1., 0.)),
        );

        let world = World::new(vec![sphere1, sphere2, lower, upper], Some(light));

        let ray = Ray::new(Tuple::point(0., 0., 0.), Tuple::vector(0., 1., 0.));

        // should terminate successfully without entering infinite recursion
        let _ = world.color_at(&ray).unwrap();
    }
}
