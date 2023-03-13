const EPSILON: f64 = 1e-5;

pub fn equal(a: f64, b: f64) -> bool {
    let c = a - b;
    c.abs() < EPSILON
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn two_floats_equal() {
        let a = 1.0;
        let b = 1.0;
        assert!(equal(a, b));
    }
}
