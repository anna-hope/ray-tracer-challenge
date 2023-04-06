use crate::{Color, Tuple};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Light {
    pub position: Tuple,
    pub intensity: Color,
}

impl Light {
    pub fn new(position: Tuple, intensity: Color) -> Self {
        Self {
            position,
            intensity,
        }
    }
}

impl Default for Light {
    fn default() -> Self {
        let position = Tuple::point(0., 0., 0.);
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
        assert_eq!(light.position, Tuple::point(0., 0., 0.));
        assert_eq!(light.intensity, Color::new(1., 1., 1.));
    }
}
