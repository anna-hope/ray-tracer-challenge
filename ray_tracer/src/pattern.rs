use crate::{color::Color, shape::Shape, Matrix, Result, Tuple};

#[derive(Debug, Clone, PartialEq)]
pub struct StripePattern {
    a: Color,
    b: Color,
    transformation: Matrix,
}

impl StripePattern {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transformation: Matrix::identity(),
        }
    }

    pub fn stripe_at(&self, point: Tuple) -> Color {
        if point.x.floor() % 2. == 0. {
            self.a
        } else {
            self.b
        }
    }

    /// Computes the stripe color for the given object at the given point.
    pub fn stripe_at_object(&self, object: &dyn Shape, world_point: Tuple) -> Result<Color> {
        // convert the world point to object space
        let object_point = object.transformation().inverse()? * world_point;

        // convert the point to pattern space
        let pattern_point = self.transformation.inverse()? * object_point;
        Ok(self.stripe_at(pattern_point))
    }

    pub fn with_transformation(mut self, transformation: Matrix) -> Self {
        self.transformation = transformation;
        self
    }
}

#[cfg(test)]
mod tests {
    use crate::{shape::sphere::Sphere, Matrix};

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

    #[test]
    fn stripes_with_object_transformation() {
        let object = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
        let pattern = StripePattern::new(Color::white(), Color::black());
        let color = pattern
            .stripe_at_object(&object, Tuple::point(1.5, 0., 0.))
            .unwrap();
        assert_eq!(color, Color::white());
    }

    #[test]
    fn stripes_with_pattern_transformation() {
        let object = Sphere::new();
        let pattern = StripePattern {
            a: Color::white(),
            b: Color::black(),
            transformation: Matrix::scaling(2., 2., 2.),
        };
        let color = pattern
            .stripe_at_object(&object, Tuple::point(1.5, 0., 0.))
            .unwrap();
        assert_eq!(color, Color::white());
    }

    #[test]
    fn stripes_with_both_object_and_pattern_transformation() {
        let object = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
        let pattern = StripePattern {
            a: Color::white(),
            b: Color::black(),
            transformation: Matrix::translation(0.5, 0., 0.),
        };

        let color = pattern
            .stripe_at_object(&object, Tuple::point(2.5, 0., 0.))
            .unwrap();
        assert_eq!(color, Color::white());
    }
}
