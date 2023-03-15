use std::ops::{Index, Mul};

use thiserror::Error;

use crate::common::equal;
use crate::Tuple;

type Result<T> = std::result::Result<T, MatrixError>;

#[derive(Debug, Error)]
pub enum MatrixError {
    #[error("Attempted to calculate an inverse of a non-invertible matrix")]
    NonInvertible,
}

#[derive(Clone, Debug)]
pub struct Matrix {
    data: Vec<Vec<f64>>,
}

impl Matrix {
    /// Does not check whether the input data has a valid shape.
    pub fn new(input_data: &[&[f64]]) -> Self {
        let data = input_data
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>();
        Self { data }
    }

    pub fn identity() -> Self {
        Self::new(&[
            &[1.0, 0.0, 0.0, 0.0],
            &[0.0, 1.0, 0.0, 0.0],
            &[0.0, 0.0, 1.0, 0.0],
            &[0.0, 0.0, 0.0, 1.0],
        ])
    }

    fn dim(&self) -> (usize, usize) {
        (self.data.len(), self.data[0].len())
    }

    fn zeros(rows: usize, columns: usize) -> Self {
        let mut data = vec![];
        for col in 0..rows {
            data.push(vec![]);
            for _ in 0..columns {
                data[col].push(0.0);
            }
        }

        Self { data }
    }

    pub fn transpose(&self) -> Self {
        // There is a much more efficient way of doing this if we only assume square matrices
        // And probably a more efficient way of doing it even with non-square matrices
        // But this works ¯\_(ツ)_/¯
        let (m, n) = self.dim();

        let mut transposed = Matrix::zeros(n, m);

        for (i, row) in self.data.iter().enumerate() {
            for (j, item) in row.iter().enumerate() {
                transposed.data[j][i] = *item;
            }
        }

        transposed
    }

    pub fn determinant(&self) -> f64 {
        if self.dim() == (2, 2) {
            let a = self[0][0];
            let b = self[0][1];
            let c = self[1][0];
            let d = self[1][1];
            a * d - b * c
        } else {
            let row = &self[0];
            let mut det = 0.0;
            for (j, val) in row.into_iter().enumerate() {
                let cofactor = self.cofactor(0, j);
                det += val * cofactor;
            }
            det
        }
    }

    pub fn submatrix(&self, row_index: usize, column_index: usize) -> Self {
        let mut data = vec![];
        for (i, row) in self.data.iter().enumerate() {
            if i != row_index {
                let mut new_row = vec![];
                for (j, val) in row.iter().enumerate() {
                    if j != column_index {
                        new_row.push(*val);
                    }
                }
                data.push(new_row);
            }
        }
        Self { data }
    }

    pub fn minor(&self, row_index: usize, column_index: usize) -> f64 {
        let submatrix = self.submatrix(row_index, column_index);
        submatrix.determinant()
    }

    pub fn cofactor(&self, row_index: usize, column_index: usize) -> f64 {
        let minor = self.minor(row_index, column_index);
        // if row + column is odd, then we negate the minor
        // otherwise, we return it as is
        if (row_index + column_index) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }

    pub fn is_invertible(&self) -> bool {
        self.determinant() != 0.0
    }

    pub fn inverse(&self) -> Result<Self> {
        if !self.is_invertible() {
            return Err(MatrixError::NonInvertible);
        }

        let (m, n) = self.dim();
        let mut new_matrix = Matrix::zeros(m, n);

        let det = self.determinant();

        for row in 0..m {
            for column in 0..n {
                let cofactor = self.cofactor(row, column);
                new_matrix.data[column][row] = cofactor / det;
            }
        }

        Ok(new_matrix)
    }
}

impl Index<usize> for Matrix {
    type Output = Vec<f64>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl PartialEq for Matrix {
    fn eq(&self, other: &Self) -> bool {
        for (this_row, other_row) in self.data.iter().zip(other.data.iter()) {
            for (this_val, other_val) in this_row.into_iter().zip(other_row) {
                if !equal(*this_val, *other_val) {
                    return false;
                }
            }
        }
        true
    }
}

impl Mul<Self> for Matrix {
    type Output = Self;

    /// Only defined for 4x4 matrices, as per the book.
    /// Will likely panic for matrices of all other dimensions.
    fn mul(self, other: Self) -> Self {
        let mut result_rows = Vec::with_capacity(4);
        for row in 0..4 {
            let mut result_row = Vec::with_capacity(4);
            for col in 0..4 {
                let result = self[row][0] * other[0][col]
                    + self[row][1] * other[1][col]
                    + self[row][2] * other[2][col]
                    + self[row][3] * other[3][col];
                result_row.push(result);
            }
            result_rows.push(result_row);
        }
        Self { data: result_rows }
    }
}

impl Mul<Tuple> for Matrix {
    type Output = Tuple;

    /// Only defined for 4x4 matrices, as per the book.
    /// Will likely panic for matrices of all other dimensions.
    fn mul(self, rhs: Tuple) -> Self::Output {
        let mut results = Vec::with_capacity(4);
        for row in 0..4 {
            let result = self[row][0] * rhs.x
                + self[row][1] * rhs.y
                + self[row][2] * rhs.z
                + self[row][3] * rhs.w;
            results.push(result);
        }
        Tuple::new(results[0], results[1], results[2], results[3])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_4x4_matrix() {
        let m = Matrix::new(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[5.5, 6.5, 7.5, 8.5],
            &[9.0, 10.0, 11.0, 12.0],
            &[13.5, 14.5, 15.5, 16.5],
        ]);

        assert_eq!(m[0][0], 1.0);
        assert_eq!(m[0][3], 4.0);
        assert_eq!(m[1][0], 5.5);
        assert_eq!(m[1][2], 7.5);
        assert_eq!(m[2][2], 11.0);
        assert_eq!(m[3][0], 13.5);
        assert_eq!(m[3][2], 15.5);
    }

    #[test]
    fn build_2x2_matrix() {
        let m = Matrix::new(&[&[-3.0, 5.0], &[1.0, -2.0]]);
        assert_eq!(m[0][0], -3.0);
        assert_eq!(m[0][1], 5.0);
        assert_eq!(m[1][0], 1.0);
        assert_eq!(m[1][1], -2.0);
    }

    #[test]
    fn build_3x3_matrix() {
        let m = Matrix::new(&[&[-3.0, 5.0, 0.0], &[1.0, -2.0, 7.0], &[0.0, 1.0, 1.0]]);
        assert_eq!(m[0][0], -3.0);
        assert_eq!(m[1][1], -2.0);
        assert_eq!(m[2][2], 1.0);
    }

    #[test]
    fn identical_matrices_are_equal() {
        let a = Matrix::new(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[5.0, 6.0, 7.0, 8.0],
            &[9.0, 8.0, 7.0, 6.0],
            &[5.0, 4.0, 3.0, 2.0],
        ]);
        let b = Matrix::new(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[5.0, 6.0, 7.0, 8.0],
            &[9.0, 8.0, 7.0, 6.0],
            &[5.0, 4.0, 3.0, 2.0],
        ]);
        assert_eq!(a, b);
    }

    #[test]
    fn different_matrices_arent_equal() {
        let a = Matrix::new(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[5.0, 6.0, 7.0, 8.0],
            &[9.0, 8.0, 7.0, 6.0],
            &[5.0, 4.0, 3.0, 2.0],
        ]);
        let b = Matrix::new(&[
            &[2.0, 3.0, 4.0, 5.0],
            &[6.0, 7.0, 8.0, 9.0],
            &[8.0, 7.0, 6.0, 5.0],
            &[4.0, 3.0, 2.0, 1.0],
        ]);
        assert_ne!(a, b);
    }

    #[test]
    fn multiply_matrices() {
        let a = Matrix::new(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[5.0, 6.0, 7.0, 8.0],
            &[9.0, 8.0, 7.0, 6.0],
            &[5.0, 4.0, 3.0, 2.0],
        ]);

        let b = Matrix::new(&[
            &[-2.0, 1.0, 2.0, 3.0],
            &[3.0, 2.0, 1.0, -1.0],
            &[4.0, 3.0, 6.0, 5.0],
            &[1.0, 2.0, 7.0, 8.0],
        ]);
        let expected = Matrix::new(&[
            &[20.0, 22.0, 50.0, 48.0],
            &[44.0, 54.0, 114.0, 108.0],
            &[40.0, 58.0, 110.0, 102.0],
            &[16.0, 26.0, 46.0, 42.0],
        ]);
        assert_eq!(a * b, expected);
    }

    #[test]
    fn multiply_matrix_by_tuple() {
        let a = Matrix::new(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[2.0, 4.0, 4.0, 2.0],
            &[8.0, 6.0, 4.0, 1.0],
            &[0.0, 0.0, 0.0, 1.0],
        ]);
        let b = Tuple::new(1.0, 2.0, 3.0, 1.0);
        let expected = Tuple::new(18.0, 24.0, 33.0, 1.0);
        assert_eq!(a * b, expected);
    }

    #[test]
    fn multiply_by_identity_matrix() {
        let a = Matrix::new(&[
            &[0.0, 1.0, 2.0, 4.0],
            &[1.0, 2.0, 4.0, 8.0],
            &[2.0, 4.0, 8.0, 16.0],
            &[4.0, 8.0, 16.0, 32.0],
        ]);
        let expected = a.clone();
        assert_eq!(a * Matrix::identity(), expected);
    }

    #[test]
    fn multiply_identity_matrix_by_tuple() {
        let a = Tuple::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(Matrix::identity() * a, a);
    }

    #[test]
    fn transpose_matrix() {
        let a = Matrix::new(&[
            &[0.0, 9.0, 3.0, 0.0],
            &[9.0, 8.0, 0.0, 8.0],
            &[1.0, 8.0, 5.0, 3.0],
            &[0.0, 0.0, 5.0, 8.0],
        ]);
        let a_transpose = Matrix::new(&[
            &[0.0, 9.0, 1.0, 0.0],
            &[9.0, 8.0, 8.0, 0.0],
            &[3.0, 0.0, 5.0, 5.0],
            &[0.0, 8.0, 3.0, 8.0],
        ]);
        assert_eq!(a.transpose(), a_transpose);
    }

    #[test]
    fn transpose_identity_matrix() {
        let a = Matrix::identity();
        assert_eq!(a.transpose(), Matrix::identity());
    }

    #[test]
    fn determinant_2x2_matrix() {
        let a = Matrix::new(&[&[1.0, 5.0], &[-3.0, 2.0]]);
        assert_eq!(a.determinant(), 17.0);
    }

    #[test]
    fn submatrix_3x3_matrix_is_2x2_matrix() {
        let a = Matrix::new(&[&[1.0, 5.0, 0.0], &[-3.0, 2.0, 7.0], &[0.0, 6.0, -3.0]]);
        let sub_a = Matrix::new(&[&[-3.0, 2.0], &[0.0, 6.0]]);
        assert_eq!(a.submatrix(0, 2), sub_a);
    }

    #[test]
    fn submatrix_4x4_matrix_is_3x3_matrix() {
        let a = Matrix::new(&[
            &[-6., 1., 1., 6.],
            &[-8., 5., 8., 6.],
            &[-1., 0., 8., 2.],
            &[-7., 1., -1., 1.],
        ]);
        let sub_a = Matrix::new(&[&[-6.0, 1.0, 6.0], &[-8.0, 8.0, 6.0], &[-7.0, -1.0, 1.0]]);
        assert_eq!(a.submatrix(2, 1), sub_a);
    }

    #[test]
    fn minor_3x3_matrix() {
        let a = Matrix::new(&[&[3.0, 5.0, 0.0], &[2.0, -1.0, -7.0], &[6.0, -1.0, 5.0]]);
        let b = a.submatrix(1, 0);
        let det_b = b.determinant();
        assert_eq!(det_b, 25.0);
        assert_eq!(a.minor(1, 0), 25.0);
    }

    #[test]
    fn cofactor_3x3_matrix() {
        let a = Matrix::new(&[&[3.0, 5.0, 0.0], &[2.0, -1.0, -7.0], &[6.0, -1.0, 5.0]]);
        assert_eq!(a.minor(0, 0), -12.0);
        assert_eq!(a.cofactor(0, 0), -12.0);
        assert_eq!(a.minor(1, 0), 25.0);
        assert_eq!(a.cofactor(1, 0), -25.0);
    }

    #[test]
    fn determinant_3x3_matrix() {
        let a = Matrix::new(&[&[1.0, 2.0, 6.0], &[-5.0, 8.0, -4.0], &[2.0, 6.0, 4.0]]);
        assert_eq!(a.cofactor(0, 0), 56.0);
        assert_eq!(a.cofactor(0, 1), 12.0);
        assert_eq!(a.determinant(), -196.0);
    }

    #[test]
    fn determinant_4x4_matrix() {
        let a = Matrix::new(&[
            &[-2., -8., 3., 5.],
            &[-3., 1., 7., 3.],
            &[1., 2., -9., 6.],
            &[-6., 7., 7., -9.],
        ]);
        assert_eq!(a.cofactor(0, 0), 690.0);
        assert_eq!(a.cofactor(0, 1), 447.0);
        assert_eq!(a.cofactor(0, 2), 210.0);
        assert_eq!(a.cofactor(0, 3), 51.0);
        assert_eq!(a.determinant(), -4071.0);
    }

    #[test]
    fn invertible_matrix_is_invertible() {
        let a = Matrix::new(&[
            &[6., 4., 4., 4.],
            &[5., 5., 7., 6.],
            &[4., -9., 3., -7.],
            &[9., 1., 7., -6.],
        ]);
        assert_eq!(a.determinant(), -2120.0);
        assert!(a.is_invertible());
    }

    #[test]
    fn non_invertible_matrix_is_not_invertible() {
        let a = Matrix::new(&[
            &[-4., 2., -2., -3.],
            &[9., 6., 2., 6.],
            &[0., -5., 1., -5.],
            &[0., 0., 0., 0.],
        ]);
        assert_eq!(a.determinant(), 0.0);
        assert!(!a.is_invertible());
    }

    #[test]
    fn calculate_inverse() {
        let a = Matrix::new(&[
            &[-5., 2., 6., -8.],
            &[1., -5., 1., 8.],
            &[7., 7., -6., -7.],
            &[1., -3., 7., 4.],
        ]);
        let b = a.inverse().unwrap();
        assert_eq!(a.determinant(), 532.0);
        assert_eq!(a.cofactor(2, 3), -160.0);
        assert!(equal(b[3][2], -160.0 / 532.0));
        assert_eq!(a.cofactor(3, 2), 105.0);
        assert!(equal(b[2][3], 105.0 / 532.0));

        let expected_b = Matrix::new(&[
            &[0.21805, 0.45113, 0.24060, -0.04511],
            &[-0.80827, -1.45677, -0.44361, 0.52068],
            &[-0.07895, -0.22368, -0.05263, 0.19737],
            &[-0.52256, -0.81391, -0.30075, 0.30639],
        ]);
        assert_eq!(b, expected_b);
    }

    #[test]
    fn calculate_another_inverse() {
        let a = Matrix::new(&[
            &[8., -5., 9., 2.],
            &[7., 5., 6., 1.],
            &[-6., 0., 9., 6.],
            &[-3., 0., -9., -4.],
        ]);
        let expected = Matrix::new(&[
            &[-0.15385, -0.15385, -0.28205, -0.53846],
            &[-0.07692, 0.12308, 0.02564, 0.03077],
            &[0.35897, 0.35897, 0.43590, 0.92308],
            &[-0.69231, -0.69231, -0.76923, -1.92308],
        ]);

        assert_eq!(a.inverse().unwrap(), expected);
    }

    #[test]
    fn calculate_yet_another_inverse() {
        let a = Matrix::new(&[
            &[9., 3., 0., 9.],
            &[-5., -2., -6., -3.],
            &[-4., 9., 6., 4.],
            &[-7., 6., 6., 2.],
        ]);
        let expected = Matrix::new(&[
            &[-0.04074, -0.07778, 0.14444, -0.22222],
            &[-0.07778, 0.03333, 0.36667, -0.33333],
            &[-0.02901, -0.14630, -0.10926, 0.12963],
            &[0.17778, 0.06667, -0.26667, 0.33333],
        ]);
        assert_eq!(a.inverse().unwrap(), expected);
    }

    #[test]
    fn multiplying_a_product_by_inverse_gives_original_matrix() {
        let a = Matrix::new(&[
            &[3., -9., 7., 3.],
            &[3., -8., 2., -9.],
            &[-4., 4., 4., 1.],
            &[-6., 5., -1., 1.],
        ]);
        let b = Matrix::new(&[
            &[8., 2., 2., 2.],
            &[3., -1., 7., 0.],
            &[7., 0., 5., 4.],
            &[6., -2., 0., 5.],
        ]);
        let c = a.clone() * b.clone();
        assert_eq!(c * b.inverse().unwrap(), a);
    }
}
