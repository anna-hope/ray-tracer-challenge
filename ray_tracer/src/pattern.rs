use crate::{color::Color, Tuple};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StripePattern {
    a: Color,
    b: Color,
}

impl StripePattern {
    pub fn new(a: Color, b: Color) -> Self {
        Self { a, b }
    }

    pub fn stripe_at(&self, point: Tuple) -> Color {
        if point.x.floor() % 2. == 0. {
            self.a
        } else {
            self.b
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_stripe_pattern() {
        let black = Color::black();
        let white = Color::white();
        let pattern = StripePattern::new(white, black);
        assert_eq!(pattern.a, white);
        assert_eq!(pattern.b, black);
    }

    #[test]
    fn stripe_pattern_is_constant_in_y() {
        let black = Color::black();
        let white = Color::white();
        let pattern = StripePattern::new(white, black);
        assert_eq!(pattern.stripe_at(Tuple::point(0., 1., 0.)), white);
        assert_eq!(pattern.stripe_at(Tuple::point(0., 2., 0.)), white);
    }

    #[test]
    fn stripe_pattern_is_constant_in_z() {
        let black = Color::black();
        let white = Color::white();
        let pattern = StripePattern::new(white, black);
        assert_eq!(pattern.stripe_at(Tuple::point(0., 0., 1.)), white);
        assert_eq!(pattern.stripe_at(Tuple::point(0., 0., 2.)), white);
    }

    #[test]
    fn stripe_pattern_alternates_in_x() {
        let black = Color::black();
        let white = Color::white();
        let pattern = StripePattern::new(white, black);
        assert_eq!(pattern.stripe_at(Tuple::point(0., 0., 0.)), white);
        assert_eq!(pattern.stripe_at(Tuple::point(0.9, 0., 0.)), white);
        assert_eq!(pattern.stripe_at(Tuple::point(1., 0., 0.)), black);
        assert_eq!(pattern.stripe_at(Tuple::point(-0.1, 0., 0.)), black);
        assert_eq!(pattern.stripe_at(Tuple::point(-1., 0., 0.)), black);
        assert_eq!(pattern.stripe_at(Tuple::point(-1.1, 0., 0.)), white);
    }
}
