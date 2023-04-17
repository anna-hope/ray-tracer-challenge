use crate::{color::Color, shape::Shape, Matrix, Point, Result};
use std::fmt;

#[derive(Debug, PartialEq)]
pub enum PatternType {
    Stripe,
    Gradient,
    Ring,
    Checker,
    Radial,
    Solid,
    Blended,
    Perturbed,
    Test,
}

pub trait PatternClone {
    fn clone_box(&self) -> Box<dyn Pattern>;
}

pub trait Pattern: PatternClone + Send + Sync {
    fn transformation(&self) -> Matrix;

    /// Converts a world point into a color for the given shape.
    fn pattern_at_shape(&self, shape: &dyn Shape, world_point: Point) -> Result<Color> {
        // convert the world point to object space
        let object_point = shape.transformation().inverse()? * world_point;

        // convert the point to pattern space
        let pattern_point = self.transformation().inverse()? * object_point;
        Ok(self.pattern_at(pattern_point))
    }

    fn pattern_at(&self, point: Point) -> Color;

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
    use super::{solid::SolidPattern, *};

    #[derive(Debug, Clone)]
    pub struct StripePattern {
        a: Box<dyn Pattern>,
        b: Box<dyn Pattern>,
        transformation: Matrix,
    }

    impl StripePattern {
        pub fn new(a: Box<dyn Pattern>, b: Box<dyn Pattern>, transformation: Matrix) -> Self {
            Self {
                a,
                b,
                transformation,
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
                a: Box::<SolidPattern>::default(),
                b: Box::new(SolidPattern::new(Color::black())),
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

        fn pattern_at(&self, point: Point) -> Color {
            let pattern = if point.x.floor() % 2. == 0. {
                &self.a
            } else {
                &self.b
            };
            pattern.pattern_at(point)
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::{shape::sphere::Sphere, Matrix};

        use super::*;

        #[test]
        fn create_stripe_pattern() {
            let _pattern = StripePattern::default();
        }

        #[test]
        fn stripe_pattern_is_constant_in_y() {
            let white = Color::white();
            let pattern = StripePattern::default();
            assert_eq!(pattern.pattern_at(Point::new(0., 1., 0.)), white);
            assert_eq!(pattern.pattern_at(Point::new(0., 2., 0.)), white);
        }

        #[test]
        fn stripe_pattern_is_constant_in_z() {
            let white = Color::white();
            let pattern = StripePattern::default();
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 1.)), white);
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 2.)), white);
        }

        #[test]
        fn stripe_pattern_alternates_in_x() {
            let black = Color::black();
            let white = Color::white();
            let pattern = StripePattern::default();
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), white);
            assert_eq!(pattern.pattern_at(Point::new(0.9, 0., 0.)), white);
            assert_eq!(pattern.pattern_at(Point::new(1., 0., 0.)), black);
            assert_eq!(pattern.pattern_at(Point::new(-0.1, 0., 0.)), black);
            assert_eq!(pattern.pattern_at(Point::new(-1., 0., 0.)), black);
            assert_eq!(pattern.pattern_at(Point::new(-1.1, 0., 0.)), white);
        }

        #[test]
        fn stripes_with_object_transformation() {
            let object = Sphere::default().with_transformation(Matrix::scaling(2., 2., 2.));
            let pattern = StripePattern::default();
            let color = pattern
                .pattern_at_shape(&object, Point::new(1.5, 0., 0.))
                .unwrap();
            assert_eq!(color, Color::white());
        }

        #[test]
        fn stripes_with_pattern_transformation() {
            let object = Sphere::default();
            let pattern = StripePattern {
                a: Box::<SolidPattern>::default(),
                b: Box::new(SolidPattern::new(Color::black())),
                transformation: Matrix::scaling(2., 2., 2.),
            };
            let color = pattern
                .pattern_at_shape(&object, Point::new(1.5, 0., 0.))
                .unwrap();
            assert_eq!(color, Color::white());
        }

        #[test]
        fn stripes_with_both_object_and_pattern_transformation() {
            let object = Sphere::default().with_transformation(Matrix::scaling(2., 2., 2.));
            let pattern =
                StripePattern::default().with_transformation(Matrix::translation(0.5, 0., 0.));

            let color = pattern
                .pattern_at_shape(&object, Point::new(2.5, 0., 0.))
                .unwrap();
            assert_eq!(color, Color::white());
        }
    }
}

pub mod gradient {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct GradientPattern {
        pub a: Color,
        pub b: Color,
        pub transformation: Matrix,
    }

    impl GradientPattern {
        pub fn new(a: Color, b: Color, transformation: Matrix) -> Self {
            Self {
                a,
                b,
                transformation,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }
    }

    impl Pattern for GradientPattern {
        fn pattern_at(&self, point: Point) -> Color {
            let distance = self.b - self.a;
            self.a + distance * point.x
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
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(
                pattern.pattern_at(Point::new(0.25, 0., 0.)),
                Color::new(0.75, 0.75, 0.75)
            );
            assert_eq!(
                pattern.pattern_at(Point::new(0.5, 0., 0.)),
                Color::new(0.5, 0.5, 0.5)
            );
            assert_eq!(
                pattern.pattern_at(Point::new(0.75, 0., 0.)),
                Color::new(0.25, 0.25, 0.25)
            );
        }
    }
}

pub mod ring {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct RingPattern {
        pub a: Color,
        pub b: Color,
        pub transformation: Matrix,
    }

    impl RingPattern {
        pub fn new(a: Color, b: Color, transformation: Matrix) -> Self {
            Self {
                a,
                b,
                transformation,
            }
        }
    }

    impl Pattern for RingPattern {
        fn pattern_at(&self, point: Point) -> Color {
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
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(1., 0., 0.)), Color::black());
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 1.)), Color::black());
            // 0.708 = just slightly more than sqrt(2) / 2
            assert_eq!(
                pattern.pattern_at(Point::new(0.708, 0., 0.708)),
                Color::black()
            );
        }
    }
}

pub mod checker {

    use super::{solid::SolidPattern, *};

    #[derive(Debug, Clone)]
    pub struct CheckerPattern {
        pub a: Box<dyn Pattern>,
        pub b: Box<dyn Pattern>,
        pub transformation: Matrix,
    }

    impl CheckerPattern {
        pub fn new(a: Box<dyn Pattern>, b: Box<dyn Pattern>, transformation: Matrix) -> Self {
            Self {
                a,
                b,
                transformation,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }
    }

    impl Default for CheckerPattern {
        fn default() -> Self {
            Self {
                a: Box::<SolidPattern>::default(),
                b: Box::new(SolidPattern::new(Color::black())),
                transformation: Matrix::identity(),
            }
        }
    }

    impl Pattern for CheckerPattern {
        fn pattern_at(&self, point: Point) -> Color {
            if (point.x.floor() + point.y.floor() + point.z.floor()) % 2. == 0. {
                self.a.pattern_at(point)
            } else {
                self.b.pattern_at(point)
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
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(0.99, 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(1.01, 0., 0.)), Color::black());
        }

        #[test]
        fn checkers_should_repeat_in_y() {
            let pattern = CheckerPattern::default();
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(0., 0.99, 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(0., 1.01, 0.)), Color::black());
        }

        #[test]
        fn checkers_should_repeat_in_z() {
            let pattern = CheckerPattern::default();
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.99)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 1.01)), Color::black());
        }
    }
}

pub mod radial_gradient {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct RadialGradientPattern {
        a: Color,
        b: Color,
        transformation: Matrix,
    }

    impl RadialGradientPattern {
        pub fn new(a: Color, b: Color, transformation: Matrix) -> Self {
            Self {
                a,
                b,
                transformation,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }
    }

    impl Default for RadialGradientPattern {
        fn default() -> Self {
            Self {
                a: Color::white(),
                b: Color::black(),
                transformation: Matrix::identity(),
            }
        }
    }

    impl Pattern for RadialGradientPattern {
        fn pattern_type(&self) -> PatternType {
            PatternType::Radial
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        /// Interpolates between two colors radially.
        fn pattern_at(&self, point: Point) -> Color {
            let distance = (point.x.powi(2) + point.z.powi(2)).sqrt();
            let fraction = distance - distance.floor();
            self.a + (self.b - self.a) * fraction
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn radial_gradient_should_linearly_interpolate_between_colors() {
            // smoothly flow between two colors in the x direction
            let pattern = RadialGradientPattern::default();
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(
                pattern.pattern_at(Point::new(0.25, 0., 0.)),
                Color::new(0.75, 0.75, 0.75)
            );
            assert_eq!(
                pattern.pattern_at(Point::new(0.5, 0., 0.)),
                Color::new(0.5, 0.5, 0.5)
            );
        }
    }
}

pub mod solid {
    use super::*;

    #[derive(Debug, Clone, Copy)]
    pub struct SolidPattern {
        color: Color,
    }

    impl SolidPattern {
        pub fn new(color: Color) -> Self {
            Self { color }
        }
    }

    impl Default for SolidPattern {
        fn default() -> Self {
            Self {
                color: Color::white(),
            }
        }
    }

    impl Pattern for SolidPattern {
        /// Always returns the identity matrix.
        fn transformation(&self) -> Matrix {
            Matrix::identity()
        }

        fn pattern_type(&self) -> PatternType {
            PatternType::Solid
        }

        /// Always returns the same color, regardless of the point.
        fn pattern_at(&self, _point: Point) -> Color {
            self.color
        }

        /// Always returns the same color, regardless of the shape or point.
        fn pattern_at_shape(&self, _shape: &dyn Shape, _world_point: Point) -> Result<Color> {
            Ok(self.color)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn solid_pattern_always_returns_same_color() {
            let pattern = SolidPattern {
                color: Color::white(),
            };
            assert_eq!(pattern.pattern_at(Point::new(0., 0., 0.)), Color::white());
            assert_eq!(pattern.pattern_at(Point::new(1., 0., 0.)), Color::white());
        }
    }
}

pub mod blended {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct BlendedPattern {
        a: Box<dyn Pattern>,
        b: Box<dyn Pattern>,
        transformation: Matrix,
    }

    impl BlendedPattern {
        pub fn new(a: Box<dyn Pattern>, b: Box<dyn Pattern>, transformation: Matrix) -> Self {
            Self {
                a,
                b,
                transformation,
            }
        }

        pub fn with_transformation(mut self, transformation: Matrix) -> Self {
            self.transformation = transformation;
            self
        }
    }

    impl Pattern for BlendedPattern {
        fn pattern_at(&self, point: Point) -> Color {
            let color1 = self.a.pattern_at(point);
            let color2 = self.b.pattern_at(point);
            let red = (color1.red + color2.red) / 2.;
            let green = (color1.green + color2.green) / 2.;
            let blue = (color1.blue + color2.blue) / 2.;
            Color::new(red, green, blue)
        }

        fn pattern_type(&self) -> PatternType {
            PatternType::Blended
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }
    }
}

pub mod perturbed {

    use super::*;
    use perlin_noise::PerlinNoise;

    #[derive(Debug, Clone)]
    pub struct PerturbedPattern {
        pattern: Box<dyn Pattern>,
        transformation: Matrix,
    }

    impl PerturbedPattern {
        pub fn new(pattern: Box<dyn Pattern>, transformation: Matrix) -> Self {
            Self {
                pattern,
                transformation,
            }
        }
    }

    impl Pattern for PerturbedPattern {
        /// Perturbs the given point with Perlin noise.
        fn pattern_at(&self, point: Point) -> Color {
            let perlin = PerlinNoise::new();
            let noise = perlin.get3d([point.x, point.y, point.z]);
            let x = noise * 0.2 + point.x;
            let y = noise * 0.2 + point.y;
            let z = noise * 0.2 + point.z;
            self.pattern.pattern_at(Point::new(x, y, z))
        }

        fn transformation(&self) -> Matrix {
            self.transformation.clone()
        }

        fn pattern_type(&self) -> PatternType {
            PatternType::Perturbed
        }
    }
}

// this has to be pub because we use the TestPattern in tests for world.rs
#[cfg(test)]
pub mod tests {

    use crate::shape::sphere::Sphere;

    use super::*;

    #[derive(Clone)]
    pub struct TestPattern {
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

        fn pattern_at(&self, point: Point) -> Color {
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
        let test_pattern: Box<dyn Pattern> = Box::<TestPattern>::default();
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
        let shape = Sphere::default().with_transformation(Matrix::scaling(2., 2., 2.));
        let pattern: Box<dyn Pattern> = Box::<TestPattern>::default();
        let color = pattern
            .pattern_at_shape(&shape, Point::new(2., 3., 4.))
            .unwrap();
        assert_eq!(color, Color::new(1., 1.5, 2.));
    }

    #[test]
    fn pattern_with_pattern_transformation() {
        let shape = Sphere::default();
        let pattern: Box<dyn Pattern> =
            Box::new(TestPattern::default().with_transformation(Matrix::scaling(2., 2., 2.)));
        let color = pattern
            .pattern_at_shape(&shape, Point::new(2., 3., 4.))
            .unwrap();
        assert_eq!(color, Color::new(1., 1.5, 2.));
    }

    #[test]
    fn pattern_with_both_object_and_pattern_transformation() {
        let shape = Sphere::default().with_transformation(Matrix::scaling(2., 2., 2.));
        let pattern: Box<dyn Pattern> =
            Box::new(TestPattern::default().with_transformation(Matrix::translation(0.5, 1., 1.5)));
        let color = pattern
            .pattern_at_shape(&shape, Point::new(2.5, 3., 3.5))
            .unwrap();
        assert_eq!(color, Color::new(0.75, 0.5, 0.25));
    }
}
