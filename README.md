# ga-core

Geometric algebra lets you do 3D math without matrices. Rotors don't gimbal lock. Reflections compose into rotations. The geometric product unifies dot products and cross products into one operation.

This crate implements Cl(3,1) — the spacetime algebra of 4D spacetime (3 space + 1 time dimension), plus conformal embedding for unified treatment of points, lines, planes, and spheres.

```rust
use ga_core::*;
```

---

## 1. A Multivector: 16 Components, 5 Grades

In Cl(3,1), every element is a multivector with 16 components organized by grade:

```
Grade 0:  1 scalar    (numbers)
Grade 1:  4 vectors   (e₀, e₁, e₂, e₃)
Grade 2:  6 bivectors (e₀₁, e₀₂, e₀₃, e₁₂, e₁₃, e₂₃)
Grade 3:  4 trivectors (oriented volumes)
Grade 4:  1 pseudoscalar (e₀₁₂₃)
```

```rust
use ga_core::Multivector;

// A scalar (grade 0)
let s = Multivector::scalar(5.0);
println!("Scalar: grade-0 norm = {:.1}", s.grade_norm(0)); // 5.0

// A vector (grade 1) in Cl(3,1): [e₀, e₁, e₂, e₃]
// Think: [time, x, y, z] or [w, x, y, z]
let v = Multivector::vector([1.0, 2.0, 3.0, 4.0]);
println!("Vector: grade-1 norm = {:.1}", v.grade_norm(1)); // 1+2+3+4 = 10
let vp = v.vector_part();
println!("Components: e₀={:.0}, e₁={:.0}, e₂={:.0}, e₃={:.0}",
    vp[0], vp[1], vp[2], vp[3]);

// A bivector (grade 2): represents a plane of rotation
let b = Multivector::bivector([1.0, 0.0, 0.0, 0.5, 0.0, 0.0]);
// Nonzero: e₀₁=1.0 (time-x plane), e₁₂=0.5 (x-y plane)
let bp = b.bivector_part();
println!("Bivector: e₀₁={:.1}, e₁₂={:.1}", bp[0], bp[3]);

// The zero multivector
let z = Multivector::zero();
println!("Zero: is_zero? {}", z.is_zero(1e-10));
```

---

## 2. The Geometric Product: Two Vectors → Scalar + Bivector

The geometric product `ab` = `a·b` (inner) + `a∧b` (outer). One operation gives you both the dot product and the cross product.

```rust
use ga_core::Multivector;

// Two spacetime vectors
let a = Multivector::vector([1.0, 2.0, 3.0, 0.0]); // space-like
let b = Multivector::vector([0.0, 1.0, 0.0, 0.0]); // pure x-direction

// Geometric product: ab = a·b + a∧b
let ab = a.geometric_product(&b);
println!("Geometric product a*b:");
println!("  Scalar part (a·b):  {:.1}", ab.scalar_part());
// e₁*e₁ = 1 (metric), so 2*1 = 2.0

println!("  Bivector part (a∧b):");
let bv = ab.bivector_part();
println!("    e₀₁={:.1}, e₀₂={:.1}, e₀₃={:.1}", bv[0], bv[1], bv[2]);
println!("    e₁₂={:.1}, e₁₃={:.1}, e₂₃={:.1}", bv[3], bv[4], bv[5]);
// Nonzero bivector components = the plane spanned by a and b

// Norm squared: |a|² = a·a
let norm_sq = a.norm_squared();
println!("|a|² = {:.1}", norm_sq);
// 1² + 2² + 3² = 14 (spacelike: e₁²=e₂²=+1)
```

### The Metric: Cl(3,1)

```rust
use ga_core::Multivector;

// In Cl(3,1): e₀²=+1, e₁²=+1, e₂²=+1, e₃²=-1
// e₃ is the "time" direction (negative signature)

let e3 = Multivector::vector([0.0, 0.0, 0.0, 1.0]);
let e3_sq = e3.geometric_product(&e3);
println!("e₃² = {:.1}", e3_sq.scalar_part()); // -1.0 (time-like!)

let e1 = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
let e1_sq = e1.geometric_product(&e1);
println!("e₁² = {:.1}", e1_sq.scalar_part()); // +1.0 (space-like)

// Orthogonal vectors: their product is purely bivector (no scalar part)
let e1v = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
let e2v = Multivector::vector([0.0, 0.0, 1.0, 0.0]);
let wedge = e1v.wedge(&e2v);
println!("e₁∧e₂ scalar part: {:.1}", wedge.scalar_part()); // 0.0
println!("e₁∧e₂ bivector e₁₂: {:.1}", wedge.bivector_part()[3]); // 1.0

// Parallel vectors: their product is purely scalar (no bivector part)
let inner = e1v.inner(&e2v);
println!("e₁·e₂ = {:.1}", inner.scalar_part()); // 0.0 (orthogonal!)
```

---

## 3. Rotors: Rotation Without Matrices, Without Gimbal Lock

A rotor is `cos(θ/2) - sin(θ/2)B` where B is the bivector of the rotation plane. Rotation happens via the sandwich product: `v' = R v R̃`.

```rust
use ga_core::Rotor;
use std::f64::consts::FRAC_PI_2;

// Rotate 90° around the z-axis
let rz = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);

// Rotate the x-axis unit vector
let v = [1.0, 0.0, 0.0];
let rotated = rz.apply(v);
println!("Rotated [1,0,0] by 90° around z: [{:.3}, {:.3}, {:.3}]",
    rotated[0], rotated[1], rotated[2]);
// ≈ [0, 1, 0] — the x-axis became the y-axis

// Rotate around a different axis
let rx = Rotor::from_axis_angle([1.0, 0.0, 0.0], FRAC_PI_2);
let v2 = [0.0, 1.0, 0.0];
let rotated2 = rx.apply(v2);
println!("Rotated [0,1,0] by 90° around x: [{:.3}, {:.3}, {:.3}]",
    rotated2[0], rotated2[1], rotated2[2]);
// ≈ [0, 0, 1] — y-axis became z-axis
```

### Compose Rotations: Multiply Rotors

```rust
use ga_core::Rotor;
use std::f64::consts::FRAC_PI_2;

// Two 90° rotations: z then x
let rz = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
let rx = Rotor::from_axis_angle([1.0, 0.0, 0.0], FRAC_PI_2);

// Compose: R = Rx * Rz (apply Rz first, then Rx)
let combined = rx.compose(&rz);

let v = [1.0, 0.0, 0.0];
let result = combined.apply(v);
println!("90° around z, then 90° around x: [{:.3}, {:.3}, {:.3}]",
    result[0], result[1], result[2]);

// Extract rotation matrix (if you need it for legacy code)
let mat = combined.to_rotation_matrix();
println!("3×3 rotation matrix:");
for row in &mat {
    println!("  [{:.3}, {:.3}, {:.3}]", row[0], row[1], row[2]);
}
```

### SLERP: Smooth Rotation Interpolation

```rust
use ga_core::Rotor;
use std::f64::consts::PI;

let start = Rotor::identity();
let end = Rotor::from_axis_angle([0.0, 0.0, 1.0], PI);

// Interpolate: t=0 is start, t=1 is end
for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
    let interp = start.slerp(&end, t);
    let rotated = interp.apply([1.0, 0.0, 0.0]);
    println!("t={:.2}: [1,0,0] → [{:.3}, {:.3}, 0]",
        t, rotated[0], rotated[1]);
}
// Smooth arc from [1,0,0] through [0,1,0] to [-1,0,0]
```

---

## 4. Reflection and Projection

```rust
use ga_core::Conformal;

let point = [3.0, 4.0, 5.0];

// Reflect through the yz-plane (normal = x-axis)
let reflected = Conformal::reflect(point, [1.0, 0.0, 0.0], 0.0);
println!("Reflected [{:.0},{:.0},{:.0}] through yz-plane: [{:.0},{:.0},{:.0}]",
    point[0], point[1], point[2],
    reflected[0], reflected[1], reflected[2]);
// [-3, 4, 5] — x component flipped

// Reflect through an arbitrary plane: z=2
let reflected2 = Conformal::reflect(point, [0.0, 0.0, 1.0], 2.0);
println!("Reflected through z=2: [{:.1},{:.1},{:.1}]",
    reflected2[0], reflected2[1], reflected2[2]);
// [3, 4, -1] — reflected about z=2: 5→2-(5-2)=-1

// Distance is preserved by reflection
let d1 = (point[0].powi(2) + point[1].powi(2) + point[2].powi(2)).sqrt();
let d2 = (reflected[0].powi(2) + reflected[1].powi(2) + reflected[2].powi(2)).sqrt();
println!("Distance preserved: {:.3} == {:.3}", d1, d2);

// Project onto a plane
let projected = Conformal::project_onto_plane(point, [0.0, 0.0, 1.0], 0.0);
println!("Project onto z=0: [{:.1},{:.1},{:.1}]",
    projected[0], projected[1], projected[2]);
// [3, 4, 0] — drop the z component
```

---

## 5. Conformal Embedding: 3D → 5D

The conformal model embeds Euclidean 3D space into 5D conformal space. This buys you: rotations become linear operations, translations become rotations, and distances become inner products.

```rust
use ga_core::{Conformal, Multivector, Rotor};

// Embed a 3D point into conformal space
let p = [3.0, 4.0, 0.0];
let embedded = Conformal::embed_point(p);

println!("Conformal embedding of ({},{},{}):", p[0], p[1], p[2]);
println!("  e₊ component: {:.1}", embedded.c[1]); // 1.0 (origin)
println!("  spatial: ({:.1}, {:.1}, {:.1})", embedded.c[2], embedded.c[3], embedded.c[4]);
println!("  ½|p|² = {:.1}", embedded.scalar_part()); // 0.5*(9+16) = 12.5

// Extract it back
let extracted = Conformal::extract_point(&embedded);
println!("Roundtrip: ({:.1}, {:.1}, {:.1}) — lossless", extracted[0], extracted[1], extracted[2]);

// Conformal distance: inner product encodes distance²
let p1 = Conformal::embed_point([1.0, 0.0, 0.0]);
let p2 = Conformal::embed_point([4.0, 0.0, 0.0]);
let dist = Conformal::conformal_distance(&p1, &p2);
println!("Distance between (1,0,0) and (4,0,0): {:.1}", dist); // 3.0

// Rotate in conformal space
let rotor = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
let rotated = Conformal::rotate([1.0, 0.0, 0.0], &rotor);
println!("Conformal rotation: [{:.3}, {:.3}, {:.3}]", rotated[0], rotated[1], rotated[2]);
// ≈ [0, 1, 0]
```

### Midpoint and Barycenter in Conformal Space

```rust
use ga_core::Conformal;

// Midpoint
let mid = Conformal::midpoint([0.0, 0.0, 0.0], [4.0, 2.0, 6.0]);
println!("Midpoint: [{:.1}, {:.1}, {:.1}]", mid[0], mid[1], mid[2]);
// [2, 1, 3]

// Barycenter with weights
let points = vec![[0.0, 0.0, 0.0], [6.0, 0.0, 0.0], [0.0, 6.0, 0.0]];
let weights = vec![1.0, 1.0, 1.0];
let centroid = Conformal::barycenter(&points, &weights);
println!("Centroid of triangle: [{:.2}, {:.2}, {:.2}]",
    centroid[0], centroid[1], centroid[2]);
// [2, 2, 0]

// Weighted barycenter (closer to first point)
let weights2 = vec![3.0, 1.0, 1.0];
let weighted = Conformal::barycenter(&points, &weights2);
println!("Weighted centroid: [{:.2}, {:.2}, {:.2}]",
    weighted[0], weighted[1], weighted[2]);
// [1.2, 1.2, 0]
```

### Plane as a Multivector

```rust
use ga_core::Conformal;

// A plane: z = 5 (parallel to xy-plane, offset by 5)
let plane = Conformal::plane([0.0, 0.0, 1.0], 5.0);
println!("Plane: normal=({:.0},{:.0},{:.0}), d={:.0}",
    plane.c[2], plane.c[3], plane.c[4], -plane.scalar_part());

// The plane itself is a multivector — you can operate on it
// with the same geometric product, reflections, etc.
```

---

## 6. The Reverse and the Conjugate

```rust
use ga_core::Multivector;

let b = Multivector::bivector([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

// Reverse: flips sign of grade-2 and grade-3 components
let rev = b.reverse();
println!("Original bivector: {:?}", b.bivector_part());
println!("Reversed bivector: {:?}", rev.bivector_part());
// [−1, −2, −3, −4, −5, −6] — all bivector signs flipped
// The reverse is used in: |M|² = M̃M

// Conjugate: negates grade-1 and grade-3
let v = Multivector::vector([1.0, 2.0, 3.0, 4.0]);
let conj = v.conjugate();
println!("Conjugate of vector: {:?}", conj.vector_part());
// [-1, -2, -3, -4]

// Dual: multiply by the pseudoscalar
let v2 = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
let dual = v2.dual();
println!("Dual of e₀:");
println!("  Pseudoscalar component: {:.1}", dual.c[15]); // maps to highest grade
// The dual converts vectors to trivectors (and vice versa)
// In 3D: dual of a vector = a plane, dual of a bivector = a normal vector
```

---

## 7. Full Example: Rigid Body Rotation

```rust
use ga_core::{Rotor, Conformal};
use std::f64::consts::PI;

// A spacecraft at position (10, 0, 5)
let position = [10.0, 0.0, 5.0];

// It has three axis vectors (its orientation frame)
let forward = [1.0, 0.0, 0.0];  // points where it's going
let up      = [0.0, 0.0, 1.0];  // points up
let right   = [0.0, 1.0, 0.0];  // points to starboard

// Yaw 45° (turn left around up-axis)
let yaw = Rotor::from_axis_angle([0.0, 0.0, 1.0], PI / 4.0);

// Pitch 30° (nose down, around right-axis)
let pitch = Rotor::from_axis_angle([0.0, 1.0, 0.0], PI / 6.0);

// Combine: yaw then pitch
let orientation = pitch.compose(&yaw);

// Apply to the frame vectors
let new_forward = orientation.apply(forward);
let new_up = orientation.apply(up);
let new_right = orientation.apply(right);

println!("After yaw 45° + pitch 30°:");
println!("  Forward: [{:.3}, {:.3}, {:.3}]", new_forward[0], new_forward[1], new_forward[2]);
println!("  Up:      [{:.3}, {:.3}, {:.3}]", new_up[0], new_up[1], new_up[2]);
println!("  Right:   [{:.3}, {:.3}, {:.3}]", new_right[0], new_right[1], new_right[2]);

// No gimbal lock. No matrix multiplications. No quaternion renormalization.
// Just compose rotors. Apply them. Done.
```

---

## Why Geometric Algebra?

| Problem | Traditional | Geometric Algebra |
|---------|------------|-------------------|
| Rotate a vector | 3×3 matrix or quaternion | `R v R̃` |
| Compose rotations | Matrix multiply or q×q | Rotor compose (same thing, cleaner) |
| Gimbal lock | Euler angles break | Rotors never lose a degree of freedom |
| Dot product | `a·b` | `½(ab + ba)` scalar part |
| Cross product | `a×b` (only in 3D!) | `a∧b` bivector (works in ANY dimension) |
| Reflect point | Explicit formula | `−nan` sandwich |
| Norm squared | `Σxᵢ²` | `v ṽ` scalar part |
| Normal to plane | Cross product of edges | Dual of the plane bivector |

**One product (geometric), one operation (sandwich), all of 3D math.**

---

## API Reference

| Type | What it does |
|------|-------------|
| `Multivector` | The fundamental element. 16 components in Cl(3,1). Geometric product, wedge, inner, reverse, dual. |
| `Rotor` | Even-grade element for rotations. From axis-angle, compose, SLERP, sandwich product, rotation matrix. |
| `Conformal` | 3D→5D embedding. Points, planes, reflection, projection, midpoint, barycenter, conformal distance. |
