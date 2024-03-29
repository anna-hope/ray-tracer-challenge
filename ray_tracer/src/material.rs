use crate::{color::Color, light::Light, pattern::Pattern, shape::Shape, Point, Result, Vector};

#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub pattern: Option<Box<dyn Pattern>>,
    pub reflectivity: f64,
    pub transparency: f64,
    pub refractive_index: f64,
    pub casts_shadow: bool,
}

impl Material {
    /// Adds together the material's ambient, diffuse, and specular components,
    /// weighted by the angles between the eye_vector and the normal_vector.
    pub fn lighting(
        &self,
        object: &dyn Shape,
        light: &dyn Light,
        point: Point,
        eye_vector: Vector,
        normal_vector: Vector,
        light_intensity: f64,
    ) -> Result<Color> {
        let color = if let Some(pattern) = &self.pattern {
            pattern.pattern_at_shape(object, point)?
        } else {
            self.color
        };

        // combine the surface color with the light's color/intensity
        let effective_color = color * light.intensity();

        // find the direction to the light surface
        let light_vector = (light.position() - point).norm();

        // compute the ambient contribution
        let ambient = effective_color * self.ambient;

        // light_dot_normal represents the cosine of the angle between
        // the light vector and the normal vector. A negative number means
        // the light is on the other side of the surface.
        let light_dot_normal = light_vector.dot(normal_vector);

        let (diffuse, specular) = if light_dot_normal < 0. {
            // set both to black because the light is obstructed
            (Color::black(), Color::black())
        } else {
            // compute the diffuse contribution
            let diffuse = effective_color * self.diffuse * light_dot_normal * light_intensity;

            // reflect_dot_eye represents the cosine of the angle between
            // the reflection vector and the eye vector. A negative number means
            // the light reflects away from the eye.
            let reflect_vector = (-light_vector).reflect(normal_vector);
            let reflect_dot_eye = reflect_vector.dot(eye_vector);

            let specular = if reflect_dot_eye <= 0. {
                // set it to black because it reflects away from the eye
                Color::black() * light_intensity
            } else {
                // compute the specular contribution
                let factor = reflect_dot_eye.powf(self.shininess);
                light.intensity() * self.specular * factor * light_intensity
            };

            (diffuse, specular)
        };

        Ok(ambient + specular + diffuse)
    }
}

impl Default for Material {
    /// Initialize material with the default parameters:
    /// color = white, ambient = 0.1, diffuse = 0.9, specular = 0.9,
    /// shininess = 200, pattern = None, reflective = 0., transparency: 0.,
    /// refractive_index: 1.
    fn default() -> Self {
        Self {
            color: Color::white(),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.,
            pattern: None,
            reflectivity: 0.,
            transparency: 0.,
            refractive_index: 1.,
            casts_shadow: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        light::PointLight, pattern::stripe::StripePattern, shape::sphere::Sphere, world::World,
    };
    use std::sync::Arc;

    use super::*;

    #[test]
    fn default_material() {
        let material = Material::default();
        assert_eq!(material.color, Color::white());
        assert_eq!(material.ambient, 0.1);
        assert_eq!(material.diffuse, 0.9);
        assert_eq!(material.specular, 0.9);
        assert_eq!(material.shininess, 200.);
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface() {
        let material = Material::default();
        let position = Point::new(0., 0., 0.);
        let eye_vector = Vector::new(0., 0., -1.);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let object = Sphere::default();
        let result = material
            .lighting(&object, &light, position, eye_vector, normal_vector, 1.)
            .unwrap();
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface_eye_offset_45() {
        let material = Material::default();
        let position = Point::new(0., 0., 0.);
        let val = 2.0_f64.sqrt() / 2.0;
        let eye_vector = Vector::new(0., val, -val);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let object = Sphere::default();
        let result = material
            .lighting(&object, &light, position, eye_vector, normal_vector, 1.)
            .unwrap();
        assert_eq!(result, Color::white());
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45() {
        let material = Material::default();
        let position = Point::new(0., 0., 0.);
        let eye_vector = Vector::new(0., 0., -1.);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::white());
        let object = Sphere::default();
        let result = material
            .lighting(&object, &light, position, eye_vector, normal_vector, 1.)
            .unwrap();
        let val = 0.7364;
        assert_eq!(result, Color::new(val, val, val));
    }

    #[test]
    fn ligthing_with_eye_in_path_of_reflection_vector() {
        let material = Material::default();
        let position = Point::new(0., 0., 0.);
        let val = 2.0_f64.sqrt() / 2.;
        let eye_vector = Vector::new(0., -val, -val);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::white());
        let object = Sphere::default();
        let result = material
            .lighting(&object, &light, position, eye_vector, normal_vector, 1.)
            .unwrap();
        let val2 = 1.6364;
        assert_eq!(result, Color::new(val2, val2, val2));
    }

    #[test]
    fn lighting_with_eye_behind_surface() {
        let material = Material::default();
        let position = Point::new(0., 0., 0.);
        let eye_vector = Vector::new(0., 0., -1.);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., 10.), Color::white());
        let object = Sphere::default();
        let result = material
            .lighting(&object, &light, position, eye_vector, normal_vector, 1.)
            .unwrap();
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn ligthing_with_surface_in_shadow() {
        let material = Material::default();
        let position = Point::new(0., 0., 0.);
        let eye_vector = Vector::new(0., 0., -1.);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let object = Sphere::default();
        let result = material
            .lighting(&object, &light, position, eye_vector, normal_vector, 0.)
            .unwrap();
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_pattern_applied() {
        let material = Material {
            ambient: 1.,
            diffuse: 0.,
            specular: 0.,
            pattern: Some(Box::<StripePattern>::default()),
            ..Default::default()
        };

        let eye_vector = Vector::new(0., 0., -1.);
        let normal_vector = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let shape = Sphere::default();
        let color1 = material
            .lighting(
                &shape,
                &light,
                Point::new(0.9, 0., 0.),
                eye_vector,
                normal_vector,
                1.,
            )
            .unwrap();

        let color2 = material
            .lighting(
                &shape,
                &light,
                Point::new(1.1, 0., 0.),
                eye_vector,
                normal_vector,
                1.,
            )
            .unwrap();

        assert_eq!(color1, Color::white());
        assert_eq!(color2, Color::black());
    }

    #[test]
    fn reflectivity_default_material() {
        let material = Material::default();
        assert_eq!(material.reflectivity, 0.);
    }

    #[test]
    fn transparency_refractive_index_for_default_material() {
        let material = Material::default();
        assert_eq!(material.transparency, 0.);
        assert_eq!(material.refractive_index, 1.);
    }

    #[test]
    fn lighting_uses_light_intensity_to_attenuate_color() {
        let mut world = World::default();
        world.lights = vec![Box::new(PointLight::new(
            Point::new(0., 0., -10.),
            Color::white(),
        ))];

        let shape = Arc::get_mut(&mut world.objects[0]).unwrap();
        let material = Material {
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.,
            color: Color::white(),
            ..Default::default()
        };
        shape.set_material(material);

        let point = Point::new(0., 0., -1.);
        let eye_vector = Vector::new(0., 0., -1.);
        let normal_vector = Vector::new(0., 0., -1.);

        let examples = [
            (1., Color::white()),
            (0.5, Color::new(0.55, 0.55, 0.55)),
            (0., Color::new(0.1, 0.1, 0.1)),
        ];

        for (light_intensity, expected) in examples {
            let result = shape
                .material()
                .lighting(
                    shape,
                    &*world.lights[0],
                    point,
                    eye_vector,
                    normal_vector,
                    light_intensity,
                )
                .unwrap();
            assert_eq!(result, expected);
        }
    }
}
