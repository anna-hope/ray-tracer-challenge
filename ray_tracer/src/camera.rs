use crate::{canvas::Canvas, world::World, Matrix, Ray, Result, Tuple};

use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;

#[derive(Debug)]
pub struct Camera {
    hsize: usize,
    vsize: usize,
    transformation: Matrix,
    half_width: f64,
    half_height: f64,
    pixel_size: f64,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, field_of_view: f64) -> Self {
        let half_view = (field_of_view / 2.).tan();
        let aspect = hsize as f64 / vsize as f64;

        let (half_width, half_height) = if aspect >= 1. {
            (half_view, half_view / aspect)
        } else {
            (half_view * aspect, half_view)
        };

        let pixel_size = half_width * 2. / hsize as f64;

        Self {
            hsize,
            vsize,
            transformation: Matrix::identity(),
            half_width,
            half_height,
            pixel_size,
        }
    }

    pub fn with_transformation(mut self, transformation: Matrix) -> Self {
        self.transformation = transformation;
        self
    }

    /// Computes the world coordinates at the center of the given pixel,
    /// then constructs a ray that passes through that point.
    fn ray_for_pixel(&self, px: usize, py: usize) -> Result<Ray> {
        // the offset from the edge of the canvas to the pixel's center
        let x_offset = (px as f64 + 0.5) * self.pixel_size;
        let y_offset = (py as f64 + 0.5) * self.pixel_size;

        // the untransformed coordinates of the pixel in the world space
        // (the camera looks toward -z, so +x is the the *left*)
        let world_x = self.half_width - x_offset;
        let world_y = self.half_height - y_offset;

        // using the camera matrix, transform the canvas point and the origin,
        // and then compute the ray's direction vector
        // (the canvas is at z=-1)
        let pixel = self.transformation.inverse()? * Tuple::point(world_x, world_y, -1.);
        let origin = self.transformation.inverse()? * Tuple::point(0., 0., 0.);
        let direction = (pixel - origin).norm();

        Ok(Ray::new(origin, direction))
    }

    /// Creates a canvas and casts a ray through each of its pixels,
    /// coloring the pixels with the colors of the corresponding intersections.
    /// Side-effect: creates and displays a progress bar to stderr.
    pub fn render(&self, world: &World) -> Result<Canvas> {
        let mut image = Canvas::new(self.hsize, self.vsize);

        let style = ProgressStyle::with_template(
            "{msg} {elapsed:>5} -- {eta:5} {bar:40.cyan/blue} {pos:>7}/{len:7} {percent}%",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar());

        let colors = (0..self.vsize - 1)
            .into_par_iter()
            .progress_with_style(style)
            .with_message("Rendering...")
            .map(|y| {
                let mut colors = vec![];
                for x in 0..self.hsize - 1 {
                    let ray = self.ray_for_pixel(x, y)?;
                    let color = world.color_at(&ray)?;
                    colors.push(color);
                }
                Ok(colors)
            })
            .collect::<Result<Vec<_>>>()?;

        for (y, row_colors) in colors.iter().enumerate() {
            for (x, color) in row_colors.iter().enumerate() {
                image.write_pixel(x, y, color.to_owned());
            }
        }
        Ok(image)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{equal, transformation::compute_view_transformation, world::World, Color};
    use std::f64::consts::PI;

    #[test]
    fn construct_camera() {
        let hsize = 160;
        let vsize = 120;
        let field_of_view = PI / 2.;
        let camera = Camera::new(hsize, vsize, field_of_view);
        assert_eq!(camera.hsize, hsize);
        assert_eq!(camera.vsize, vsize);
        assert_eq!(camera.transformation, Matrix::identity());
    }

    #[test]
    fn pixel_size_horizontal_canvas() {
        let camera = Camera::new(200, 125, PI / 2.);
        assert!(equal(camera.pixel_size, 0.01));
    }

    #[test]
    fn pixel_size_vertical_canvas() {
        let camera = Camera::new(125, 200, PI / 2.);
        assert!(equal(camera.pixel_size, 0.01));
    }

    #[test]
    fn construct_ray_through_center_of_canvas() {
        let camera = Camera::new(201, 101, PI / 2.);
        let ray = camera.ray_for_pixel(100, 50).unwrap();
        assert_eq!(ray.origin, Tuple::point(0., 0., 0.));
        assert_eq!(ray.direction, Tuple::vector(0., 0., -1.));
    }

    #[test]
    fn construct_ray_through_corner_of_canvas() {
        let camera = Camera::new(201, 101, PI / 2.);
        let ray = camera.ray_for_pixel(0, 0).unwrap();
        assert_eq!(ray.origin, Tuple::point(0., 0., 0.));
        assert_eq!(ray.direction, Tuple::vector(0.66519, 0.33259, -0.66851));
    }

    #[test]
    fn construct_ray_when_camera_is_transformed() {
        let camera = Camera::new(201, 101, PI / 2.)
            .with_transformation(Matrix::identity().translate(0., -2., 5.).rotate_y(PI / 4.));
        let ray = camera.ray_for_pixel(100, 50).unwrap();
        assert_eq!(ray.origin, Tuple::point(0., 2., -5.));

        let val = 2.0_f64.sqrt() / 2.;
        assert_eq!(ray.direction, Tuple::vector(val, 0., -val));
    }

    #[test]
    fn rendering_world_with_camera() {
        let world = World::default();
        let from = Tuple::point(0., 0., -5.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        let camera = Camera::new(11, 11, PI / 2.)
            .with_transformation(compute_view_transformation(from, to, up).unwrap());
        let image = camera.render(&world).unwrap();
        assert_eq!(image.pixel_at(5, 5), Color::new(0.38066, 0.47583, 0.2855));
    }
}
