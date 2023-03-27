use crate::{color::Color, shape::Shape, Matrix, Result, Tuple};
use std::fmt;

#[derive(Debug, PartialEq)]
pub enum PatternType {
    Stripe,
    Gradient,
    Ring,
    Checker,
    Test,
}

pub trait PatternClone {
    fn clone_box(&self) -> Box<dyn Pattern>;
}

pub trait Pattern: PatternClone {
    fn transformation(&self) -> Matrix;

    /// Converts a world point into a color for the given shape.
    fn pattern_at_shape(&self, shape: &dyn Shape, world_point: Tuple) -> Result<Color> {
        // convert the world point to object space
        let object_point = shape.transformation().inverse()? * world_point;

        // convert the point to pattern space
        let pattern_point = self.transformation().inverse()? * object_point;
        Ok(self.pattern_at(pattern_point))
    }

    fn pattern_at(&self, point: Tuple) -> Color;

    fn pattern_type(&self) -> PatternType;
}

impl fmt::Debug for dyn Pattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Pattern")
            .field("type", &self.pattern_type())
            .field("transformation", &self.transformation())
            .finish()
    }
}

impl PartialEq for dyn Pattern {
    fn eq(&self, other: &Self) -> bool {
        self.pattern_type() == other.pattern_type()
            && self.transformation() == other.transformation()
    }
}

impl<T> PatternClone for T
where
    T: 'static + Pattern + Clone,
{
    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Pattern> {
    fn clone(&self) -> Box<dyn Pattern> {
        self.clone_box()
    }
}

pub mod stripe {
    use super::*;

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

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }
    }

    impl Default for StripePattern {
        fn default() -> Self {
            Self {
                a: Color::white(),
                b: Color::black(),
                transformation: Matrix::identity(),
            }
        }
    }

    impl Pattern for StripePattern {
        fn pattern_type(&self) -> PatternType {
            PatternType::Stripe
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn pattern_at(&self, point: Tuple) -> Color {
            if point.x.floor() % 2. == 0. {
                self.a
            } else {
                self.b
            }
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
            assert_eq!(pattern.pattern_at(Tuple::point(0., 1., 0.)), white);
            assert_eq!(pattern.pattern_at(Tuple::point(0., 2., 0.)), white);
        }

        #[test]
        fn stripe_pattern_is_constant_in_z() {
            let black = Color::black();
            let white = Color::white();
            let pattern = StripePattern::new(white, black);
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 1.)), white);
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 2.)), white);
        }

        #[test]
        fn stripe_pattern_alternates_in_x() {
            let black = Color::black();
            let white = Color::white();
            let pattern = StripePattern::new(white, black);
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 0.)), white);
            assert_eq!(pattern.pattern_at(Tuple::point(0.9, 0., 0.)), white);
            assert_eq!(pattern.pattern_at(Tuple::point(1., 0., 0.)), black);
            assert_eq!(pattern.pattern_at(Tuple::point(-0.1, 0., 0.)), black);
            assert_eq!(pattern.pattern_at(Tuple::point(-1., 0., 0.)), black);
            assert_eq!(pattern.pattern_at(Tuple::point(-1.1, 0., 0.)), white);
        }

        #[test]
        fn stripes_with_object_transformation() {
            let object = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
            let pattern = StripePattern::new(Color::white(), Color::black());
            let color = pattern
                .pattern_at_shape(&object, Tuple::point(1.5, 0., 0.))
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
                .pattern_at_shape(&object, Tuple::point(1.5, 0., 0.))
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
                .pattern_at_shape(&object, Tuple::point(2.5, 0., 0.))
                .unwrap();
            assert_eq!(color, Color::white());
        }
    }
}

pub mod gradient {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct GradientPattern {
        a: Color,
        b: Color,
        transformation: Matrix,
    }

    impl Pattern for GradientPattern {
        fn pattern_at(&self, point: Tuple) -> Color {
            let distance = self.b - self.a;
            let fraction = point.x - point.x.floor();
            self.a + distance * fraction
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn pattern_type(&self) -> PatternType {
            PatternType::Gradient
        }
    }

    impl Default for GradientPattern {
        fn default() -> Self {
            Self {
                a: Color::white(),
                b: Color::black(),
                transformation: Matrix::identity(),
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn gradient_linearly_interpolates_between_colors() {
            let pattern = GradientPattern::default();
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 0.)), Color::white());
            assert_eq!(
                pattern.pattern_at(Tuple::point(0.25, 0., 0.)),
                Color::new(0.75, 0.75, 0.75)
            );
            assert_eq!(
                pattern.pattern_at(Tuple::point(0.5, 0., 0.)),
                Color::new(0.5, 0.5, 0.5)
            );
            assert_eq!(
                pattern.pattern_at(Tuple::point(0.75, 0., 0.)),
                Color::new(0.25, 0.25, 0.25)
            );
        }
    }
}

pub mod ring {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct RingPattern {
        a: Color,
        b: Color,
        transformation: Matrix,
    }

    impl Pattern for RingPattern {
        fn pattern_at(&self, point: Tuple) -> Color {
            if (point.x.powi(2) + point.z.powi(2)).sqrt().floor() % 2. == 0. {
                self.a
            } else {
                self.b
            }
        }

        fn pattern_type(&self) -> PatternType {
            PatternType::Ring
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }
    }

    impl Default for RingPattern {
        fn default() -> Self {
            Self {
                a: Color::white(),
                b: Color::black(),
                transformation: Matrix::identity(),
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn ring_should_extend_in_both_x_and_z() {
            let pattern = RingPattern::default();
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Tuple::point(1., 0., 0.)), Color::black());
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 1.)), Color::black());
            // 0.708 = just slightly more than sqrt(2) / 2
            assert_eq!(
                pattern.pattern_at(Tuple::point(0.708, 0., 0.708)),
                Color::black()
            );
        }
    }
}

pub mod checker {

    use super::*;

    #[derive(Debug, Clone)]
    pub struct CheckerPattern {
        a: Color,
        b: Color,
        transformation: Matrix,
    }

    impl Default for CheckerPattern {
        fn default() -> Self {
            Self {
                a: Color::white(),
                b: Color::black(),
                transformation: Matrix::identity(),
            }
        }
    }

    impl Pattern for CheckerPattern {
        fn pattern_at(&self, point: Tuple) -> Color {
            if point.x.floor() + point.y.floor() + point.z.floor() % 2. == 0. {
                self.a
            } else {
                self.b
            }
        }

        fn pattern_type(&self) -> PatternType {
            PatternType::Checker
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn checkers_should_repeat_in_x() {
            let pattern = CheckerPattern::default();
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 0.)), Color::white());
            assert_eq!(
                pattern.pattern_at(Tuple::point(0.99, 0., 0.)),
                Color::white()
            );
            assert_eq!(
                pattern.pattern_at(Tuple::point(1.01, 0., 0.)),
                Color::black()
            );
        }

        #[test]
        fn checkers_should_repeat_in_y() {
            let pattern = CheckerPattern::default();
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 0.)), Color::white());
            assert_eq!(
                pattern.pattern_at(Tuple::point(0., 0.99, 0.)),
                Color::white()
            );
            assert_eq!(
                pattern.pattern_at(Tuple::point(0., 1.01, 0.)),
                Color::black()
            );
        }

        #[test]
        fn checkers_should_repeat_in_z() {
            let pattern = CheckerPattern::default();
            assert_eq!(pattern.pattern_at(Tuple::point(0., 0., 0.)), Color::white());
            assert_eq!(
                pattern.pattern_at(Tuple::point(0., 0., 0.99)),
                Color::white()
            );
            assert_eq!(
                pattern.pattern_at(Tuple::point(0., 0., 1.01)),
                Color::black()
            );
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::shape::sphere::Sphere;

    use super::*;

    #[derive(Clone)]
    struct TestPattern {
        transformation: Matrix,
    }

    impl TestPattern {
        fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }
    }

    impl Pattern for TestPattern {
        fn pattern_type(&self) -> PatternType {
            PatternType::Test
        }

        fn pattern_at(&self, point: Tuple) -> Color {
            Color::new(point.x, point.y, point.z)
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }
    }

    impl Default for TestPattern {
        fn default() -> Self {
            Self {
                transformation: Matrix::identity(),
            }
        }
    }

    #[test]
    fn test_pattern() {
        // erase the type in order to test the trait object
        let test_pattern: Box<dyn Pattern> = Box::new(TestPattern::default());
        assert_eq!(test_pattern.transformation(), Matrix::identity());
    }

    #[test]
    fn assigning_pattern_transformation() {
        let test_pattern: Box<dyn Pattern> =
            Box::new(TestPattern::default().with_transformation(Matrix::translation(1., 2., 3.)));
        assert_eq!(
            test_pattern.transformation(),
            Matrix::translation(1., 2., 3.)
        );
    }

    #[test]
    fn pattern_with_object_transformation() {
        let shape = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
        let pattern: Box<dyn Pattern> = Box::new(TestPattern::default());
        let color = pattern
            .pattern_at_shape(&shape, Tuple::point(2., 3., 4.))
            .unwrap();
        assert_eq!(color, Color::new(1., 1.5, 2.));
    }

    #[test]
    fn pattern_with_pattern_transformation() {
        let shape = Sphere::new();
        let pattern: Box<dyn Pattern> =
            Box::new(TestPattern::default().with_transformation(Matrix::scaling(2., 2., 2.)));
        let color = pattern
            .pattern_at_shape(&shape, Tuple::point(2., 3., 4.))
            .unwrap();
        assert_eq!(color, Color::new(1., 1.5, 2.));
    }

    #[test]
    fn pattern_with_both_object_and_pattern_transformation() {
        let shape = Sphere::new().with_transformation(Matrix::scaling(2., 2., 2.));
        let pattern: Box<dyn Pattern> =
            Box::new(TestPattern::default().with_transformation(Matrix::translation(0.5, 1., 1.5)));
        let color = pattern
            .pattern_at_shape(&shape, Tuple::point(2.5, 3., 3.5))
            .unwrap();
        assert_eq!(color, Color::new(0.75, 0.5, 0.25));
    }
}
