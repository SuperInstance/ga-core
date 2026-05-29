# ga-core

**Conformal geometric algebra in Cl(3,1) — multivectors, rotors, conformal embeddings, and sandwich products in pure Rust.**

Implements the 16-component multivector algebra of Cl(3,1) spacetime. Supports the full geometric product (combining inner and outer products), rotors for rotations via sandwich products, and conformal model operations that embed Euclidean 3D into 5D conformal space for unified treatment of points, lines, planes, circles, and spheres.

## What This Gives You

- **16-component multivectors** — scalar + 4 vectors + 6 bivectors + 4 trivectors + pseudoscalar
- **Geometric product** — the fundamental product of GA (inner ⊕ outer)
- **Rotors** — axis-angle → rotor, compose, normalize, apply sandwich product RvR̃
- **Conformal model** — embed/extract Euclidean 3D points, reflect through planes
- **Zero dependencies** — pure Rust, no linear algebra crates

## Quick Start

```rust
use ga_core::{Multivector, Rotor, Conformal};

// Create a rotor (30° rotation around Z axis)
let rotor = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_6);

// Rotate a point
let point = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
let rotated = rotor.rotate(&point);

// Conformal embedding
let embedded = Conformal::embed_point([1.0, 2.0, 3.0]);
let reflected = Conformal::reflect([1.0, 2.0, 3.0], [0.0, 0.0, 1.0], 0.0);
```

## API Reference

| Type | Description |
|------|-------------|
| `Multivector` | 16-component Cl(3,1) element with full arithmetic |
| `Rotor` | Even-grade element for rotations, from axis-angle or multivector |
| `Conformal` | Static methods for conformal space operations |

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
ga-core = { git = "https://github.com/SuperInstance/ga-core" }
```

## How It Fits

Part of the SuperInstance ecosystem:

- **[gpu-ga-kernel](https://github.com/SuperInstance/gpu-ga-kernel)** — CUDA version for GPU throughput
- **ga-core** — Rust CPU version (this repo)

## License

MIT
