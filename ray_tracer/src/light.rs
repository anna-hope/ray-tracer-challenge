use crate::{Color, Point};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Light {
    pub position: Point,
    pub intensity: Color,
}

impl Light {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            position,
            intensity,
        }
    }
}

impl Default for Light {
    fn default() -> Self {
        let position = Point::new(0., 0., 0.);
        let intensity = Color::white();
        Self::new(position, intensity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_light_has_position_and_intensity() {
        let light = Light::default();
        assert_eq!(light.position, Point::new(0., 0., 0.));
        assert_eq!(light.intensity, Color::new(1., 1., 1.));
    }
}
