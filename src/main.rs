use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial,Polynomial};
pub mod field;
use field::{F};
use ark_ff::{Fp512, MontBackend, BitIteratorBE, One, BigInt, Fp};
pub mod csidh_fri;
use csidh_fri::{fold, round_commit};

use ark_crypto_primitives::CRHScheme;
use ark_sponge::poseidon::PoseidonConfig;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
//use ark_crypto_primitives::*;
pub mod merkle;
use merkle::{FieldMT, poseidon_parameters};
fn main() {
    // 2 + x + 4x^2 + 6x^3 + 2x^4 + 111x^5 + 75x^6 + 5x^7 + 56x^8
    let poly = DensePolynomial {
        coeffs: vec![
            F::from(2),
            F::from(1),
            F::from(4),
            F::from(6),
            F::from(2),
            F::from(111),
            F::from(75),
            F::from(5),
            F::from(56)
        ]
    };
    let modulus = <F as PrimeField>::MODULUS;
    //let eval_at_three = Polynomial::evaluate(&poly, &F::from(3));
    //let three = F::from(3);
    //let five = F::from(5);
    //let two = F::from(2);
    //assert_eq!(eval_at_three, five);
    //println!("Hello, world!, {:?}", F::from(1u8));
//let folded = fold(poly,3u8, F::from(2u8));
//let eval_check = Polynomial::evaluate(&folded, &F::from(3));
//assert_eq!(eval_check, F::from(6287));
//const MODULUS: BigInt<8> = BigInt([1982068743014369403, 14011292126959937589, 5865710692925656869, 12081687501529634055, 6556111612370143693, 12983042349969476674, 18197551657619704906, 7328639240417282495]);
//let fp: Fp<8> = Fp(MODULUS, 8);
//type ff = Fp<BigInt<8>, 8>;//
//round_commit(poly, 3, F::from(2), F::from(10), F::from(23), 16);

let params = poseidon_parameters();
let crh_a = poseidon::CRH::<F>::evaluate(&params, poly.coeffs.clone()).unwrap();

println!("Hello, world!, {:?}", crh_a);
}
