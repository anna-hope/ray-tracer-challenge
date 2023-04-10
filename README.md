# The Ray Tracer Challenge in Rust (WIP)

This reposity is a work-in-progress implementation of [The Ray Tracer Challenge](http://raytracerchallenge.com/)
in Rust.

## Build and run

This project targets Rust 1.68.0 and above. [rustup](https://rustup.rs/) is recommended for installing the Rust toolchain.

To run:

Release mode is highly recommended, as the renderer is *much* slower in debug mode.

```
cargo run --release samples/<sample>.yaml
```

## Run tests

```
cargo test
```

## Implementation notes

I am keeping some notes about my implementation in **NOTES.md**. Please check there for more details.