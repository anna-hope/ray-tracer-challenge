# The Ray Tracer Challenge
## Implementation notes

### General

- f64 used for all floats
- I use `transformation` in place of `transform` (from the book) for transformation matrices
  because I think that's clearer

### Points and vectors

I implemented the functions given in the book as `point` and `vector` as associated functions on `Tuple`.
So, to get a point, you do `Tuple::point(...)`, and likewise for a vector, you do `Tuple::vector(...)`.

### Matrices

- Matrix multiplication is defined only for 4x4 matrices (and hardcoded as such), as per the book
- Matrix transformations are implemented with the help of associated functions on the `Matrix` struct

### Spheres

- I implemented object on `Intersection` as a dynamic trait object, so that I can store any object that implements the `Shape` trait as the object of the Intersection.
This is uglier than I would like, and prevents simple equality comparisons
between a `Sphere` that was downcast to a `Shape` and a `Sphere` that was not.
I tried to implement `PartialEq` for `Shape` and the concrete `Sphere` type,
but it didn't work due to lifetime issues. I'm not sure how to fix this.

- For a `Sphere` with a non-default (i.e. non-identity matrix) transformation, the book has you
instantiate a new `Sphere`, and then set the transformation on that `Sphere`. I thought this was
kind of awkward, so I added a new `with_transformation` associated function to the `Sphere` struct, which instantiates a new Sphere with the given transformation. Saves a line of code,
and no longer requires a Sphere to be mutable.

### Planes

#### Refactoring to reuse code between shapes

This is so far my least favorite part of the project (p. 117-122). The book, while not written
for a specific programming language, is by its own admission structured in an object-oriented way.
Rust is not an object-oriented language. This limits the amount of code reuse that can be done,
and what can be done feels hacky.

For example, to reuse the code in `normal_at`, we have to introduce the `local_normal_at` method
to the `Shape` trait, which becomes available to all users of structs that implement the trait.
This is not ideal, because it exposes an implementation detail to the user. The user should not
have to know that the `local_normal_at` method exists, and should not be able to call it.

With the `intersect` and `local_intersect` methods suggested in the book, it is even trickier,
because the I implemented `intersect` as part of a dedicated `Intersect` trait. That way,
we can implement `Intersect` for structs that don't implement `Shape`, like `World`.
This is a good example of composition over inheritance.
However, that makes it impossible to reuse the code in `intersect` as the book suggests, because the `Intersect` trait is independent of the `Shape` trait, and the `local_intersect` method relies on
the `transformation` method, which is part of the `Shape` trait (or equivalent field in the concrete struct).

### Patterns

This is so far my least favorite part of the codebase. Boxed dynamic trait objects
`(Box<dyn Pattern>)` galore! Working with patterns (e.g. nesting them) is unwieldy and
unergonomic --
you have to box them at every step. The whole things feels very not Rusty.

In addition, it feels like some of the patterns are not implemented correctly. For example,
the `BlendedPattern` does not seem to do much.

## Ideas for future improvements/enhancements

- The ray tracer technically supports multiple light sources, but rendering scenes with more than one light source is very slow
and seems somewhat broken.
- The YAML format to describe scenes as specified by the book doesn't include support for material patterns. I need to
  figure out how to extend it so that all the pattern code doesn't just go to waste.
