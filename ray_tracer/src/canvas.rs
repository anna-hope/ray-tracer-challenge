use crate::color::Color;

const MAX_PIXEL_VALUE: u8 = 255;
const MAX_LINE_LEN: usize = 70;

/// We use row-major order for the Vec of pixels.
pub struct Canvas {
    pub width: usize,
    pub height: usize,
    pixels: Vec<Vec<Color>>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Self {
        let mut pixels = Vec::with_capacity(height);
        for _ in 0..height {
            let mut pixels_row = Vec::with_capacity(width);
            for _ in 0..width {
                pixels_row.push(Color::default());
            }
            pixels.push(pixels_row)
        }

        Self {
            width,
            height,
            pixels,
        }
    }

    /// Panics if x or y are outside the canvas.
    pub fn pixel_at(&self, x: usize, y: usize) -> Color {
        self.pixels[y][x]
    }

    /// Panics if x or y are outside the canvas.
    pub fn write_pixel(&mut self, x: usize, y: usize, color: Color) {
        self.pixels[y][x] = color;
    }

    /// Get a PPM-format string for the pixel values of the canvas.
    pub fn to_ppm(&self) -> String {
        let mut lines = vec![];
        lines.push("P3".to_string());
        lines.push(format!("{} {}", self.width, self.height));
        lines.push(MAX_PIXEL_VALUE.to_string());

        for pixels_row in self.pixels.iter() {
            let mut line = vec![];
            for pixel in pixels_row {
                let red = scale_and_clamp(pixel.red);
                let green = scale_and_clamp(pixel.green);
                let blue = scale_and_clamp(pixel.blue);
                line.push(red);
                line.push(green);
                line.push(blue);
            }

            let mut line_strings = vec![];
            let mut num_chars_line = 0;
            for val in line {
                if num_chars_line + 4 > MAX_LINE_LEN {
                    line_strings.push("\n".to_string());
                    num_chars_line = 0;
                }
                if !line_strings.is_empty() && line_strings.last() != Some(&"\n".to_string()) {
                    line_strings.push(" ".to_string());
                    num_chars_line += 1;
                }
                line_strings.push(val.to_string());

                // 3 chars for each subpixel
                num_chars_line += 3;
            }

            lines.push(line_strings.join(""));
        }

        lines.join("\n") + "\n"
    }
}

/// Scales to 0-255, and clamps to be between those values.
fn scale_and_clamp(val: f64) -> u8 {
    let val_scaled = (val * MAX_PIXEL_VALUE as f64).round() as u8;
    val_scaled.clamp(0, 255)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_canvas() {
        let c = Canvas::new(10, 20);
        assert_eq!(c.width, 10);
        assert_eq!(c.height, 20);
        for pixels_row in c.pixels {
            for pixel in pixels_row {
                assert_eq!(pixel, Color::default());
            }
        }
    }

    #[test]
    fn write_pixels_to_canvas() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1.0, 0.0, 0.0);
        c.write_pixel(2, 3, red);
        assert_eq!(c.pixel_at(2, 3), red);
    }

    #[test]
    fn ppm_header() {
        let c = Canvas::new(5, 3);
        let ppm = c.to_ppm();
        let first_three_lines = ppm.lines().take(3).collect::<Vec<_>>();
        assert_eq!(first_three_lines[0], "P3");
        assert_eq!(first_three_lines[1], "5 3");
        assert_eq!(first_three_lines[2], "255");
    }

    #[test]
    fn ppm_pixels() {
        let mut c = Canvas::new(5, 3);
        let c1 = Color::new(1.5, 0.0, 0.0);
        let c2 = Color::new(0.0, 0.5, 0.0);
        let c3 = Color::new(-0.5, 0.0, 1.0);

        c.write_pixel(0, 0, c1);
        c.write_pixel(2, 1, c2);
        c.write_pixel(4, 2, c3);
        let ppm = c.to_ppm();
        let lines = ppm.lines().collect::<Vec<_>>();
        assert_eq!(lines[3], "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines[4], "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines[5], "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
    }

    #[test]
    fn ppm_split_long_lines() {
        let mut c = Canvas::new(10, 2);
        // iterate over every i,j index in the pixels
        // and set them all to the specified color
        for i in 0..c.pixels.len() {
            for j in 0..c.pixels[i].len() {
                c.pixels[i][j] = Color::new(1.0, 0.8, 0.6);
            }
        }

        let ppm = c.to_ppm();
        let lines = ppm.lines().collect::<Vec<_>>();
        assert_eq!(
            lines[3],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[4],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
        assert_eq!(
            lines[5],
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines[6],
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
    }

    #[test]
    fn ppm_ends_in_newline() {
      let c = Canvas::new(5, 3);
      let ppm = c.to_ppm();
      let last_char = ppm.chars().last().unwrap();
      assert_eq!(last_char, '\n');
    }
}
