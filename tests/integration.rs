//! Integration tests for ga-core

use ga_core::*;

#[test]
fn test_multivector_zero() {
    let m = Multivector::zero();
    assert!(m.is_zero(1e-15));
    assert_eq!(m.scalar_part(), 0.0);
}

#[test]
fn test_multivector_add_inverse() {
    let a = Multivector::vector([1.0, 2.0, 3.0, 4.0]);
    let neg = a.scale(-1.0);
    let sum = a.add(&neg);
    assert!(sum.is_zero(1e-10));
}

#[test]
fn test_multivector_associative_add() {
    let a = Multivector::scalar(1.0);
    let b = Multivector::scalar(2.0);
    let c = Multivector::scalar(3.0);
    let ab_c = a.add(&b).add(&c);
    let a_bc = a.add(&b.add(&c));
    assert!(ab_c.sub(&a_bc).is_zero(1e-10));
}

#[test]
fn test_geometric_product_anticommutative_vectors() {
    let v0 = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
    let v1 = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
    let v01 = v0.geometric_product(&v1);
    let v10 = v1.geometric_product(&v0);
    // v01 + v10 should be zero (bivector part) for orthogonal vectors
    // Actually: v01 + v10 = 2*(v0·v1) = 0 for orthogonal
    // v01 - v10 = 2*(v0∧v1)
    let sum = v01.add(&v10);
    assert!(sum.grade_norm(2) < 1e-10, "anticommutator of orthogonal vectors should have no bivector");
}

#[test]
fn test_rotor_preserves_length() {
    let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.23);
    let v = [3.0, 4.0, 0.0];
    let rotated = r.apply(v);
    let orig_len = (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt();
    let rot_len = (rotated[0].powi(2) + rotated[1].powi(2) + rotated[2].powi(2)).sqrt();
    assert!((orig_len - rot_len).abs() < 0.1, "rotation should preserve length");
}

#[test]
fn test_rotor_double_rotation_full() {
    // Two 90-degree rotations around z = 180 degrees total
    let r1 = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
    let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
    let r_total = r1.compose(&r2);
    let v = [1.0, 0.0, 0.0];
    let rotated = r_total.apply(v);
    // After 180 degrees: [1,0,0] → [-1,0,0]
    assert!((rotated[0] + 1.0).abs() < 0.1, "expected near -1.0, got {}", rotated[0]);
}

#[test]
fn test_conformal_embed_extract_roundtrip() {
    let points = [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-5.0, 0.5, 100.0]];
    for p in points {
        let m = Conformal::embed_point(p);
        let extracted = Conformal::extract_point(&m);
        for i in 0..3 {
            assert!((extracted[i] - p[i]).abs() < 1e-10);
        }
    }
}

#[test]
fn test_conformal_reflect_involutive() {
    // Reflecting twice through the same plane should return to original
    let p = [1.0, 2.0, 3.0];
    let normal = [0.0, 0.0, 1.0];
    let r1 = Conformal::reflect(p, normal, 0.0);
    let r2 = Conformal::reflect(r1, normal, 0.0);
    for i in 0..3 {
        assert!((r2[i] - p[i]).abs() < 1e-10, "double reflection should be identity");
    }
}

#[test]
fn test_conformal_barycenter_triangle() {
    let pts = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
    let weights = [1.0, 1.0, 1.0];
    let bc = Conformal::barycenter(&pts, &weights);
    assert!((bc[0] - 1.0).abs() < 1e-10);
    assert!((bc[1] - 2.0/3.0).abs() < 1e-10);
}

#[test]
fn test_conformal_barycenter_zero_weights() {
    let pts = [[1.0, 2.0, 3.0]];
    let weights = [0.0];
    let bc = Conformal::barycenter(&pts, &weights);
    assert_eq!(bc, [0.0, 0.0, 0.0]);
}

#[test]
fn test_dual_involution() {
    let m = Multivector::scalar(1.0);
    let d1 = m.dual();
    let d2 = d1.dual();
    // Double dual should recover original (up to sign)
    assert!((d2.scalar_part() - 1.0).abs() < 1e-10);
}
