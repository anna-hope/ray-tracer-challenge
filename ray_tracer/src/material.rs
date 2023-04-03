use crate::{color::Color, light::PointLight, pattern::Pattern, shape::Shape, Result, Tuple};

#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub pattern: Option<Box<dyn Pattern>>,
    pub reflective: f64,
    pub transparency: f64,
    pub refractive_index: f64,
}

impl Material {
    /// Adds together the material's ambient, diffuse, and specular components,
    /// weighted by the angles between the eye_vector and the normal_vector.
    pub fn lighting(
        &self,
        object: Box<dyn Shape>,
        light: PointLight,
        point: Tuple,
        eye_vector: Tuple,
        normal_vector: Tuple,
        in_shadow: bool,
    ) -> Result<Color> {
        let color = if let Some(pattern) = &self.pattern {
            pattern.pattern_at_shape(object, point)?
        } else {
            self.color
        };

        // combine the surface color with the light's color/intensity
        let effective_color = color * light.intensity;

        // find the direction to the light surface
        let light_vector = (light.position - point).norm();

        // compute the ambient contribution
        let ambient = effective_color * self.ambient;

        // light_dot_normal represents the cosine of the angle between
        // the light vector and the normal vector. A negative number means
        // the light is on the other side of the surface.
        let light_dot_normal = light_vector.dot(&normal_vector)?;

        let diffuse: Color;
        let specular: Color;

        if light_dot_normal < 0. {
            // set both to black because the light is obstructed
            diffuse = Color::default();
            specular = Color::default();
        } else {
            // compute the diffuse contribution
            diffuse = effective_color * self.diffuse * light_dot_normal;

            // reflect_dot_eye represents the cosine of the angle between
            // the reflection vector and the eye vector. A negative number means
            // the light reflects away from the eye.
            let reflect_vector = (-light_vector).reflect(normal_vector)?;
            let reflect_dot_eye = reflect_vector.dot(&eye_vector)?;

            if reflect_dot_eye <= 0. {
                // set it to black because it reflects away from the eye
                specular = Color::default();
            } else {
                // compute the specular contribution
                let factor = reflect_dot_eye.powf(self.shininess);
                specular = light.intensity * self.specular * factor;
            }
        }

        if in_shadow {
            // only the ambient light illuminates the material if we're in shadow
            Ok(ambient)
        } else {
            // add the three contributions to get the final shading
            Ok(ambient + diffuse + specular)
        }
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
            reflective: 0.,
            transparency: 0.,
            refractive_index: 1.,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{pattern::stripe::StripePattern, shape::sphere::Sphere};

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
        let position = Tuple::point(0., 0., 0.);
        let eye_vector = Tuple::vector(0., 0., -1.);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 0., -10.), Color::white());
        let object = Box::<Sphere>::default();
        let result = material
            .lighting(object, light, position, eye_vector, normal_vector, false)
            .unwrap();
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface_eye_offset_45() {
        let material = Material::default();
        let position = Tuple::point(0., 0., 0.);
        let val = 2.0_f64.sqrt() / 2.0;
        let eye_vector = Tuple::vector(0., val, -val);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 0., -10.), Color::white());
        let object = Box::<Sphere>::default();
        let result = material
            .lighting(object, light, position, eye_vector, normal_vector, false)
            .unwrap();
        assert_eq!(result, Color::white());
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45() {
        let material = Material::default();
        let position = Tuple::point(0., 0., 0.);
        let eye_vector = Tuple::vector(0., 0., -1.);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 10., -10.), Color::white());
        let object = Box::<Sphere>::default();
        let result = material
            .lighting(object, light, position, eye_vector, normal_vector, false)
            .unwrap();
        let val = 0.7364;
        assert_eq!(result, Color::new(val, val, val));
    }

    #[test]
    fn ligthing_with_eye_in_path_of_reflection_vector() {
        let material = Material::default();
        let position = Tuple::point(0., 0., 0.);
        let val = 2.0_f64.sqrt() / 2.;
        let eye_vector = Tuple::vector(0., -val, -val);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 10., -10.), Color::white());
        let object = Box::<Sphere>::default();
        let result = material
            .lighting(object, light, position, eye_vector, normal_vector, false)
            .unwrap();
        let val2 = 1.6364;
        assert_eq!(result, Color::new(val2, val2, val2));
    }

    #[test]
    fn lighting_with_eye_behind_surface() {
        let material = Material::default();
        let position = Tuple::point(0., 0., 0.);
        let eye_vector = Tuple::vector(0., 0., -1.);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 0., 10.), Color::white());
        let object = Box::<Sphere>::default();
        let result = material
            .lighting(object, light, position, eye_vector, normal_vector, false)
            .unwrap();
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn ligthing_with_surface_in_shadow() {
        let material = Material::default();
        let position = Tuple::point(0., 0., 0.);
        let eye_vector = Tuple::vector(0., 0., -1.);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 0., -10.), Color::white());
        let in_shadow = true;
        let object = Box::<Sphere>::default();
        let result = material
            .lighting(
                object,
                light,
                position,
                eye_vector,
                normal_vector,
                in_shadow,
            )
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

        let eye_vector = Tuple::vector(0., 0., -1.);
        let normal_vector = Tuple::vector(0., 0., -1.);
        let light = PointLight::new(Tuple::point(0., 0., -10.), Color::white());
        let shape = Sphere::default();
        let color1 = material
            .lighting(
                Box::new(shape.clone()),
                light,
                Tuple::point(0.9, 0., 0.),
                eye_vector,
                normal_vector,
                false,
            )
            .unwrap();

        let color2 = material
            .lighting(
                Box::new(shape),
                light,
                Tuple::point(1.1, 0., 0.),
                eye_vector,
                normal_vector,
                false,
            )
            .unwrap();

        assert_eq!(color1, Color::white());
        assert_eq!(color2, Color::black());
    }

    #[test]
    fn reflectivity_default_material() {
        let material = Material::default();
        assert_eq!(material.reflective, 0.);
    }

    #[test]
    fn transparency_refractive_index_for_default_material() {
        let material = Material::default();
        assert_eq!(material.transparency, 0.);
        assert_eq!(material.refractive_index, 1.);
    }
}
