# ga-core

Conformal geometric algebra using Cl(3,1) spacetime algebra — multivectors, rotors, and conformal embeddings.

## Usage

```rust
use ga_core::{Multivector, Rotor, Conformal};

// Create a multivector in Cl(3,1)
let v = Multivector::vector([1.0, 2.0, 3.0, 0.0]);
let scalar = Multivector::scalar(5.0);

// Geometric product
let product = v.geometric_product(&scalar);

// Rotation via rotor (sandwich product)
let rotor = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
let rotated = rotor.apply([1.0, 0.0, 0.0]);

// Conformal embedding
let embedded = Conformal::embed_point([1.0, 2.0, 3.0]);
```

## Features

- **16-component multivectors** for Cl(3,1) spacetime algebra
- **Geometric, wedge, and inner products**
- **Rotors** with axis-angle, composition, slerp, and rotation matrix extraction
- **Conformal embedding** of Euclidean 3D into 5D conformal space
- **Reflection and projection** operations

## Tests

32 tests, all passing. `cargo test` to run.

## License

MIT

Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem.
