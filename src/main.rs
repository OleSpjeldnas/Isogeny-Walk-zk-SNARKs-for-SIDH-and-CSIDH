use ark_crypto_primitives::crh::sha256::digest::generic_array::typenum::array;
use ark_ff::{MontFp};
use ark_poly::{univariate::DensePolynomial,Polynomial};
pub mod field;
use ark_std::iterable::Iterable;
use field::{FqConfig, F as Fp, Fq2 as F};
use ark_ff::{MontBackend, BitIteratorBE, One, BigInt};
pub mod matrix;
use matrix::*;
pub mod csidh_fri;
use csidh_fri::*;
use ark_ff::Field;
use ndarray::prelude::Array2;
use ndarray::*;
use ndarray_linalg::*;
use std::hash::{Hash, Hasher};
use nalgebra::*;
use ark_crypto_primitives::CRHScheme;
use ark_sponge::poseidon::PoseidonConfig;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
//use ark_crypto_primitives::*;
pub mod merkle;
use merkle::{FieldMT, poseidon_parameters, FieldPath};

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
    //let modulus = <F as PrimeField>::MODULUS;
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
let l_vec: Vec<usize> = vec![2,3,2];
let r: F = F::from(3);
let s: F = F::new(MontFp!("46125640644791984503063093131192152308860971046482148971221392549658928357566369654107368811182716603910699644514503519777239072768"), MontFp!("46125640644791984503063093131192152308860971046482148971221392549658928357566369654107368811182716603910699644514503519777239072767"));
let s_ord: u8 = 24;
//let s: F = F::new(Fp::from(17), Fp::from(0)).sqrt().unwrap();
//println!("s: {:?}", s);

//let mtrees = commit(poly, l_vec, s, r, s_ord);
//let mut mat: Array2<F> = Array::zeros((0, 10));
//let relevant_primes: Vec<u64> = vec![3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587];
//let relevant_primes: Vec<u64> = vec![4, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587];
//let tt: Vec<u64> = vec![2,5,6];
let gen = raise_to_power(s, 2u64, 210);
let gen2 = raise_to_power(gen, 3u64, 132);
let gen3 = raise_to_power(gen2, 1223u64, 1);
//let (mt, f1, )
let f1 = fold(&poly, 2, F::from(2));
let (mt1, r1, eval1) = round_commit(&poly, &gen3, &r, &s_ord);
let (mt2, r2, eval2) = round_commit(&f1, &gen3.pow(&[2]), &r.pow(&[2]), &(s_ord/2));
let index: usize = 5;
//println!("length points1:{:?}", eval1.len());
//println!("length points2:{:?}", eval2.len());
let mut points: Vec<F> = Vec::new();
points.push(eval2[5]);
points.push(eval1[5]);
points.push(eval1[17]);
let t: F = gen3.pow(&[12]);
assert_eq!(points[0], f1.evaluate(&(r.pow(&[2])*gen3.pow(&[10]))));
assert_eq!(points[1], poly.evaluate(&(r*gen3.pow(&[5]))));
assert_eq!(points[2], poly.evaluate(&(r*gen3.pow(&[5])*t)));
//println!("Correct: {:?}", verify_fold_at_index(points, gen3, t, 2, F::from(1)));
let m1: F = t/(t-F::from(1));
let m2: F = F::from(1)/(t-F::from(1));
let m3: F = F::from(1)/(gen3*(t-F::from(1)));
let m4: F = -F::from(1)/(gen3*(t-F::from(1)));

let row1: Vec<F> = vec![F::from(1), r*gen3.pow(&[5])];
let row2: Vec<F> = vec![F::from(1), r*gen3.pow(&[5])*t];
let mat: Vec<Vec<F>> = vec![row1, row2];
let y: Vec<F> = solve_linear_system(mat, vec![points[1], points[2]]);
//let (lower, upper) = lu_decomposition(&mat, 2);
//let claimed_lower = vec![vec![F::from(1), F::from(0)], vec![F::from(1), F::from(1)]];

//let sol1: F = m1*points[1] + m2*points[2];
//let sol2: F = m3*points[1] + m4*points[2];
//let c1 = lower[0][0]*upper[0][0] + lower[0][1]*upper[1][0];
//let c2 = lower[0][0]*upper[0][1] + lower[0][1]*upper[1][1];
//let c3 = lower[1][0]*upper[0][0] + lower[1][1]*upper[1][0];
//let c4 = lower[1][0]*upper[0][1] + lower[1][1]*upper[1][1];
//let mm = vec![vec![c1, c2], vec![c3, c4]];
//assert_eq!(vec![vec![m1, m2], vec![m3, m4]], inv.unwrap());
assert_eq!(y[0]+F::from(2)*y[1], points[0]);

//return;
let (paths, points, roots) = prove(poly, l_vec.clone(), gen3, r, 24, 3);
let b = verify(paths, points, roots, l_vec, gen3, r, 24, 3);


println!("Hello, world!, {:?}",b);

}
use std::collections::hash_map::DefaultHasher;
fn calculate_hash<T: Hash>(t: &T, n: u64) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish() % n
}

fn raise_to_power(x: F, v: u64, n: u8) -> F {
    let mut s = x;
    for i in 0..n{
    s = s.pow(&[v]);}
s
}