use super::*;

fn prove() {
let witness: DensePolynomial<F> = DensePolynomial::from_coefficients_vec(vec![F::from(1u128), F::from(2u128), F::from(3u128)]);
let mut coeffs: Vec<F> = Vec::new();
coe
for i in 1..witness.coeffs.len() {
    coeffs.push(witness.coeffs[i]+witness.coeffs[i-1]);
}