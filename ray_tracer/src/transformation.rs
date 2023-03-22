use crate::{Matrix, Result, Tuple};

pub fn compute_view_transformation(from: Tuple, to: Tuple, up: Tuple) -> Result<Matrix> {
    let forward = (to - from).norm();
    let up_norm = up.norm();
    let left = forward.cross(&up_norm)?;
    let true_up = left.cross(&forward)?;
    let orientation = Matrix::new(&[
        &[left.x, left.y, left.z, 0.],
        &[true_up.x, true_up.y, true_up.z, 0.],
        &[-forward.x, -forward.y, -forward.z, 0.],
        &[0., 0., 0., 1.],
    ]);
    Ok(orientation * Matrix::translation(-from.x, -from.y, -from.z))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transformation_matrix_default_orientation() {
        // the default orientation is the matrix you get
        // if your view parameters don't require any transformations
        // i.e. it's the identity matrix
        // this test has orientation looking from the origin
        // along the `z` axis in the negative direction
        // with `up` in the positive `y` direction
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., -1.);
        let up = Tuple::vector(0., 1., 0.);
        let transformation = compute_view_transformation(from, to, up).unwrap();
        assert_eq!(transformation, Matrix::identity());
    }

    #[test]
    fn view_transformation_matrix_looking_in_positive_z_direction() {
        // looking in the positive `z` direction is like
        // looking in the mirror, so the transformation is a reflection
        // which is the same as scaling everything by -1
        let from = Tuple::point(0., 0., 0.);
        let to = Tuple::point(0., 0., 1.);
        let up = Tuple::vector(0., 1., 0.);
        let transformation = compute_view_transformation(from, to, up).unwrap();
        assert_eq!(transformation, Matrix::scaling(-1., 1., -1.));
    }

    #[test]
    fn view_transformation_moves_world() {
        // the eye is 8 units along the `z` axis
        // and points back at the origin
        // this actually moves the world (and not the eye)
        // backwards 8 points along the `z` axis
        let from = Tuple::point(0., 0., 8.);
        let to = Tuple::point(0., 0., 0.);
        let up = Tuple::vector(0., 1., 0.);
        let transformation = compute_view_transformation(from, to, up).unwrap();
        assert_eq!(transformation, Matrix::translation(0., 0., -8.));
    }

    #[test]
    fn arbitrary_view_transformation() {
        // looking in some arbitrary direction
        // produces a matrix that is a combination of shearing,
        // scaling, and translation
        let from = Tuple::point(1., 3., 2.);
        let to = Tuple::point(4., -2., 8.);
        let up = Tuple::vector(1., 1., 0.);
        let transformation = compute_view_transformation(from, to, up).unwrap();
        let expected = Matrix::new(&[
            &[-0.50709, 0.50709, 0.67612, -2.36643],
            &[0.76772, 0.60609, 0.12122, -2.82843],
            &[-0.35857, 0.59761, -0.71714, 0.],
            &[0., 0., 0., 1.],
        ]);
        assert_eq!(transformation, expected);
    }
}
