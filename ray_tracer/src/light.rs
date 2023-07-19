use crate::{world::World, Color, Point, Result};

pub trait Light {
    fn intensity_at(&self, point: Point, world: &World) -> Result<f64>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            position,
            intensity,
        }
    }
}

impl Default for PointLight {
    fn default() -> Self {
        let position = Point::new(0., 0., 0.);
        let intensity = Color::white();
        Self::new(position, intensity)
    }
}

impl Light for PointLight {
    fn intensity_at(&self, point: Point, world: &World) -> Result<f64> {
        if world.is_shadowed(point, self.position)? {
            Ok(0.)
        } else {
            Ok(1.)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::equal;

    use super::*;

    #[test]
    fn point_light_has_position_and_intensity() {
        let light = PointLight::default();
        assert_eq!(light.position, Point::new(0., 0., 0.));
        assert_eq!(light.intensity, Color::new(1., 1., 1.));
    }

    #[test]
    fn point_lights_evaluate_light_intensity_at_given_point() {
        let world = World::default();
        let light = world.lights[0];

        let examples = [
            (Point::new(0., 1.0001, 0.), 1.0),
            (Point::new(-1.0001, 0., 0.), 1.0),
            (Point::new(0., 0., -1.0001), 1.0),
            (Point::new(0., 0., 1.0001), 0.),
            (Point::new(1.0001, 0., 0.), 0.),
            (Point::new(0., -1.0001, 0.), 0.),
            (Point::new(0., 0., 0.), 0.),
        ];

        for (point, expected) in examples {
            let intensity = light.intensity_at(point, &world).unwrap();
            assert!(equal(intensity, expected));
        }
    }
}
