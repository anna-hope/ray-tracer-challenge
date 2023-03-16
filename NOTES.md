# The Ray Tracer Challenge
## Implementation notes

### General

- f64 used for all floats

### Points and vectors

I implemented the functions given in the book as `point` and `vector` as associated functions on `Tuple`.
So, to get a point, you do `Tuple::point(...)`, and likewise for a vector, you do `Tuple::vector(...)`.

### Matrices

- Matrix multiplication is defined only for 4x4 matrices (and hardcoded as such), as per the book
- Matrix transformations are implemented with the help of associated functions on the `Matrix` struct
