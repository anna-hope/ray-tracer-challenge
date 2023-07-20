use crate::{world::World, Color, Point, Result, Vector};

pub trait Light: Send + Sync + 'static {
    fn intensity_at(&self, point: Point, world: &World) -> Result<f64>;

    fn position(&self) -> Point;

    fn intensity(&self) -> Color;
}

impl std::fmt::Debug for dyn Light {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Light")
            .field("position", &self.position())
            .field("intensity", &self.intensity())
            .finish()
    }
}

impl PartialEq for dyn Light {
    fn eq(&self, other: &Self) -> bool {
        self.position() == other.position() && self.intensity() == other.intensity()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointLight {
    position: Point,
    intensity: Color,
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

    fn position(&self) -> Point {
        self.position
    }

    fn intensity(&self) -> Color {
        self.intensity
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AreaLight {
    corner: Point,
    uvec: Vector,
    usteps: usize,
    vvec: Vector,
    vsteps: usize,
    intensity: Color,
    samples: usize,
    position: Point,
}

impl AreaLight {
    pub fn new(
        corner: Point,
        full_uvec: Vector,
        usteps: usize,
        full_vvec: Vector,
        vsteps: usize,
        intensity: Color,
    ) -> Self {
        let uvec = full_uvec / (usteps as f64);
        let vvec = full_vvec / (vsteps as f64);
        let samples = usteps * vsteps;
        let center = (full_uvec + full_vvec) / 2.;
        let position = Point::new(center.x, center.y, center.z);

        Self {
            corner,
            uvec,
            usteps,
            vvec,
            vsteps,
            intensity,
            samples,
            position,
        }
    }

    fn point_on_light(&self, u: usize, v: usize) -> Point {
        self.corner + self.uvec * (u as f64 + 0.5) + self.vvec * (v as f64 + 0.5)
    }
}

impl Light for AreaLight {
    fn intensity_at(&self, point: Point, world: &World) -> Result<f64> {
        let mut total = 0.;

        for v in 0..self.vsteps {
            for u in 0..self.usteps {
                let light_position = self.point_on_light(u, v);
                if !world.is_shadowed(point, light_position)? {
                    total += 1.;
                }
            }
        }

        Ok(total / (self.samples as f64))
    }

    fn position(&self) -> Point {
        self.position
    }

    fn intensity(&self) -> Color {
        self.intensity
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
        let light = &world.lights[0];

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

    #[test]
    fn create_area_light() {
        let corner = Point::new(0., 0., 0.);
        let vector1 = Vector::new(2., 0., 0.);
        let vector2 = Vector::new(0., 0., 1.);
        let light = AreaLight::new(corner, vector1, 4, vector2, 2, Color::white());

        assert_eq!(light.corner, corner);
        assert_eq!(light.uvec, Vector::new(0.5, 0., 0.));
        assert_eq!(light.usteps, 4);
        assert_eq!(light.vvec, Vector::new(0., 0., 0.5));
        assert_eq!(light.vsteps, 2);
        assert_eq!(light.samples, 8);
        assert_eq!(light.position, Point::new(1., 0., 0.5));
    }

    #[test]
    fn finding_single_point_on_area_light() {
        let corner = Point::new(0., 0., 0.);
        let vector1 = Vector::new(2., 0., 0.);
        let vector2 = Vector::new(0., 0., 1.);
        let light = AreaLight::new(corner, vector1, 4, vector2, 2, Color::white());

        let examples = [
            (0, 0, Point::new(0.25, 0., 0.25)),
            (1, 0, Point::new(0.75, 0., 0.25)),
            (0, 1, Point::new(0.25, 0., 0.75)),
            (2, 0, Point::new(1.25, 0., 0.25)),
            (3, 1, Point::new(1.75, 0., 0.75)),
        ];

        for (u, v, expected) in examples {
            let point = light.point_on_light(u, v);
            assert_eq!(point, expected);
        }
    }

    #[test]
    fn area_light_intensity_method() {
        let world = World::default();
        let corner = Point::new(-0.5, -0.5, -5.);
        let vector1 = Vector::new(1., 0., 0.);
        let vector2 = Vector::new(0., 1., 0.);
        let light = AreaLight::new(corner, vector1, 2, vector2, 2, Color::white());

        let examples = [
            (Point::new(0., 0., 2.), 0.),
            (Point::new(1., -1., 2.), 0.25),
            (Point::new(1.5, 0., 2.), 0.5),
            (Point::new(1.25, 1.25, 3.), 0.75),
            (Point::new(0., 0., -2.), 1.),
        ];

        for (point, expected) in examples {
            let intensity = light.intensity_at(point, &world).unwrap();
            assert_eq!(intensity, expected);
        }
    }
}
