//! # GA Core Tutorial
//!
//! A progressive introduction to conformal geometric algebra with Cl(3,1).
//!
//! Run with: `cargo run --example tutorial`

use ga_core::{Conformal, Multivector, Rotor};

fn main() {
    lesson_1_scalar_and_zero();
    lesson_2_vectors_and_grade();
    lesson_3_geometric_product();
    lesson_4_wedge_and_inner();
    lesson_5_rotations_with_rotors();
    lesson_6_rotor_composition_and_slerp();
    lesson_7_conformal_embedding();
    lesson_8_reflections_and_projections();
}

// ── Lesson 1: Scalars and the Zero Multivector ──────────────────────────
fn lesson_1_scalar_and_zero() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 1: Scalars and the Zero Multivector");
    println!("═══════════════════════════════════════════");
    println!();

    // The zero multivector — every component is 0.
    let z = Multivector::zero();
    println!("Zero multivector is_zero(1e-10)? {}", z.is_zero(1e-10));

    // A pure scalar has only the grade-0 component set.
    let s = Multivector::scalar(7.5);
    println!("Scalar(7.5) → scalar_part = {}", s.scalar_part());
    println!("  grade_norm(0) = {}", s.grade_norm(0));
    println!("  grade_norm(1) = {} (no vector part)", s.grade_norm(1));
    println!();

    // Basic arithmetic on scalars.
    let a = Multivector::scalar(3.0);
    let b = Multivector::scalar(4.0);
    let sum = a.add(&b);
    let doubled = a.scale(2.0);
    println!("3 + 4 = {}", sum.scalar_part());
    println!("3 × 2 = {}", doubled.scalar_part());
    println!();
}

// ── Lesson 2: Vectors and Grade Structure ────────────────────────────────
fn lesson_2_vectors_and_grade() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 2: Vectors and Grade Structure");
    println!("═══════════════════════════════════════════");
    println!();

    // Cl(3,1) vectors have 4 components: e0, e1, e2, e3.
    // Metric: e0²=e1²=e2²=+1, e3²=−1 (spacetime signature).
    let v = Multivector::vector([1.0, 2.0, 3.0, 0.0]);
    let vp = v.vector_part();
    println!("v = {}e0 + {}e1 + {}e2 + {}e3", vp[0], vp[1], vp[2], vp[3]);
    println!("  grade_norm(0) = {} (no scalar)", v.grade_norm(0));
    println!("  grade_norm(1) = {} (sum of |components|)", v.grade_norm(1));
    println!();

    // Timelike vector (negative norm² from e3).
    let timelike = Multivector::vector([0.0, 0.0, 0.0, 1.0]);
    println!("Timelike vector e3: norm² = {} (should be -1)", timelike.norm_squared());
    println!();

    // Spacelike vector.
    let spacelike = Multivector::vector([3.0, 4.0, 0.0, 0.0]);
    println!("Spacelike (3,4,0,0): norm² = {} (should be 25)", spacelike.norm_squared());
    println!();
}

// ── Lesson 3: The Geometric Product ─────────────────────────────────────
fn lesson_3_geometric_product() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 3: The Geometric Product");
    println!("═══════════════════════════════════════════");
    println!();

    // The geometric product is the fundamental operation of GA.
    // For vectors: ab = a·b + a∧b (inner + outer product).

    // Same vector squared → scalar (the metric norm).
    let e0 = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
    let e0_sq = e0.geometric_product(&e0);
    println!("e0·e0 = {} (e0² = +1)", e0_sq.scalar_part());

    let e3 = Multivector::vector([0.0, 0.0, 0.0, 1.0]);
    let e3_sq = e3.geometric_product(&e3);
    println!("e3·e3 = {} (e3² = −1)", e3_sq.scalar_part());
    println!();

    // Orthogonal vectors produce a bivector.
    let e1 = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
    let e2 = Multivector::vector([0.0, 0.0, 1.0, 0.0]);
    let e1e2 = e1.geometric_product(&e2);
    println!("e1·e2 has scalar part: {}", e1e2.scalar_part());
    println!("e1·e2 bivector part: {:?}", e1e2.bivector_part());
    println!("  (e1∧e2 = e12 bivector at index 8)");
    println!();
}

// ── Lesson 4: Wedge and Inner Products ──────────────────────────────────
fn lesson_4_wedge_and_inner() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 4: Wedge (Outer) and Inner Products");
    println!("═══════════════════════════════════════════");
    println!();

    let a = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
    let b = Multivector::vector([1.0, 1.0, 0.0, 0.0]);

    // Inner product: the symmetric part of the geometric product.
    let dot = a.inner(&b);
    println!("a · b = {} (dot product in metric)", dot.scalar_part());

    // Wedge product: the antisymmetric part — produces a bivector.
    let wedge = a.wedge(&b);
    println!("a ∧ b bivector part: {:?}", wedge.bivector_part());
    println!("  (only e01 component survives)");
    println!();

    // Parallel vectors: wedge = 0, inner ≠ 0.
    let c = Multivector::vector([2.0, 0.0, 0.0, 0.0]);
    let wedge_parallel = a.wedge(&c);
    let dot_parallel = a.inner(&c);
    println!("Parallel: a∧(2a) is_zero? {}", wedge_parallel.is_zero(1e-10));
    println!("Parallel: a·(2a) = {}", dot_parallel.scalar_part());
    println!();
}

// ── Lesson 5: Rotations with Rotors ─────────────────────────────────────
fn lesson_5_rotations_with_rotors() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 5: Rotations with Rotors");
    println!("═══════════════════════════════════════════");
    println!();

    // The identity rotor leaves vectors unchanged.
    let id = Rotor::identity();
    println!("Identity rotor is_identity? {}", id.is_identity(1e-10));
    println!("Identity applied to (1,2,3): {:?}", id.apply([1.0, 2.0, 3.0]));
    println!();

    // Rotate 90° around the z-axis.
    use std::f64::consts::FRAC_PI_2;
    let rz = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    let v = [1.0, 0.0, 0.0];
    let rotated = rz.apply(v);
    println!("Rotate (1,0,0) by 90° around z:");
    println!("  result = ({:.4}, {:.4}, {:.4})", rotated[0], rotated[1], rotated[2]);
    println!();

    // Extract rotation matrix.
    let mat = rz.to_rotation_matrix();
    println!("Rotation matrix:");
    for row in &mat {
        println!("  [{:.4}, {:.4}, {:.4}]", row[0], row[1], row[2]);
    }
    println!();
}

// ── Lesson 6: Rotor Composition and SLERP ───────────────────────────────
fn lesson_6_rotor_composition_and_slerp() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 6: Rotor Composition and SLERP");
    println!("═══════════════════════════════════════════");
    println!();

    // Compose two 90° rotations to get 180°.
    use std::f64::consts::FRAC_PI_2;
    let rz = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    let rz2 = rz.compose(&rz);
    let v = [1.0, 0.0, 0.0];
    let rotated_180 = rz2.apply(v);
    println!("Compose two 90° z-rotations → apply to (1,0,0):");
    println!("  ({:.4}, {:.4}, {:.4})", rotated_180[0], rotated_180[1], rotated_180[2]);
    println!();

    // Spherical linear interpolation (SLERP).
    let r0 = Rotor::identity();
    let r1 = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    println!("SLERP between identity and 90° z-rotation:");
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let interp = r0.slerp(&r1, t);
        let v = interp.apply([1.0, 0.0, 0.0]);
        println!("  t={:.2}: rotated (1,0,0) → ({:.4}, {:.4}, {:.4})", t, v[0], v[1], v[2]);
    }
    println!();
}

// ── Lesson 7: Conformal Embedding ────────────────────────────────────────
fn lesson_7_conformal_embedding() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 7: Conformal Embedding");
    println!("═══════════════════════════════════════════");
    println!();

    // Embed Euclidean points into conformal space.
    let p = [3.0, 4.0, 0.0];
    let cp = Conformal::embed_point(p);
    println!("Embed ({}, {}, {}):", p[0], p[1], p[2]);
    println!("  scalar_part (½|p|²) = {}", cp.scalar_part());
    println!("  vector_part = {:?}", cp.vector_part());
    println!();

    // Round-trip: extract the Euclidean point back.
    let recovered = Conformal::extract_point(&cp);
    println!("Extract: ({}, {}, {})", recovered[0], recovered[1], recovered[2]);
    println!();

    // Midpoint and barycenter.
    let mid = Conformal::midpoint([0.0, 0.0, 0.0], [4.0, 6.0, 2.0]);
    println!("Midpoint of (0,0,0) and (4,6,2): ({}, {}, {})", mid[0], mid[1], mid[2]);

    let bc = Conformal::barycenter(
        &[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        &[1.0, 1.0, 1.0],
    );
    println!("Barycenter (equal weights): ({:.4}, {:.4}, {:.4})", bc[0], bc[1], bc[2]);
    println!();
}

// ── Lesson 8: Reflections and Projections ────────────────────────────────
fn lesson_8_reflections_and_projections() {
    println!("═══════════════════════════════════════════");
    println!("Lesson 8: Reflections and Projections");
    println!("═══════════════════════════════════════════");
    println!();

    // Reflect a point through a plane.
    let point = [1.0, 2.0, 3.0];
    let reflected = Conformal::reflect(point, [1.0, 0.0, 0.0], 0.0);
    println!("Reflect (1,2,3) through yz-plane: ({}, {}, {})", reflected[0], reflected[1], reflected[2]);
    println!();

    // Project a point onto a plane.
    let projected = Conformal::project_onto_plane([1.0, 1.0, 5.0], [0.0, 0.0, 1.0], 0.0);
    println!("Project (1,1,5) onto xy-plane: ({}, {}, {})", projected[0], projected[1], projected[2]);
    println!();

    // Create a plane multivector and inspect it.
    let plane = Conformal::plane([0.0, 0.0, 1.0], 3.0);
    println!("Plane z=3: scalar_part = {}, vector = {:?}", plane.scalar_part(), plane.vector_part());
    println!();

    // Rotate with a conformal operation.
    let rotor = Rotor::from_axis_angle([0.0, 1.0, 0.0], std::f64::consts::PI);
    let rotated = Conformal::rotate([1.0, 0.0, 0.0], &rotor);
    println!("Conformal rotate (1,0,0) by 180° around y: ({:.4}, {:.4}, {:.4})",
        rotated[0], rotated[1], rotated[2]);
    println!();
}
