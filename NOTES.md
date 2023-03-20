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

### Spheres

- I implemented object on Intersection as a dynamic trait object, so that I can store any object that implements the `SceneObject` trait as the object of the Intersection. 
This is uglier than I would like, and prevents simple equality comparisons
between a Sphere that was downcast to a SceneObject and a Sphere that was not. I tried to implement PartialEq
for SceneObject and the concrete Sphere type, but it didn't work due to lifetime issues. I'm not sure how to fix this.