use ark_ff::{MontFp};
use ark_poly::{univariate::DensePolynomial,Polynomial};
//pub mod field;
pub mod test_field;
use ark_std::iterable::Iterable;
//use field::{FqConfig, F as Fp, Fq2 as F};
use test_field::{FqConfig, F as Fp, Fq2 as F};
use ark_ff::{MontBackend, BitIteratorBE, One, BigInt};
pub mod matrix;
use matrix::*;
pub mod csidh_fri;
pub mod get_roots;
pub mod isogeny_prove;
use isogeny_prove::{prove, verify};
use csidh_fri::*;
use ark_ff::Field;
use std::hash::{Hash, Hasher};
use std::ops::Div;
use ark_crypto_primitives::CRHScheme;
use ark_sponge::poseidon::PoseidonConfig;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
use std::fs::File;
use ark_std::test_rng;
use std::io::{stdin, Write};
//use ark_crypto_primitives::*;
pub mod merkle;
use merkle::{FieldMT, poseidon_parameters, FieldPath};
use get_roots::*;
use std::time::{Duration, Instant};
use crate::isogeny_prove::*;
fn main() {
    // l_vec contains all the folding factors
    let l_list: Vec<usize> = vec![vec![3; 4], vec![4; 2], vec![2]].concat();
    let mut s: F = MULT_GEN;
    
    
    for _ in 0..210 {
        s = s.pow(&[2]);
    }
    for _ in 0..131 {
        s = s.pow(&[3]);
    }
    let mut g: F = s.clone();
    for _ in 0..6 {
        g = g.pow(&[2]);
    }
    for _ in 0..2 {
        s = s.pow(&[3]);
    }
    let r: F = F::new(Fp::from(5), Fp::from(3));
    //let order_group: u64 = 2u64.pow(6);
    //let w = FFT_GEN.pow(&[order_group]).inverse().unwrap();
    // Witness polynomial
    let witness: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("new_coeffs.txt").unwrap() };
    // Witness(x+1)
    let n = witness.coeffs.len();

    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }
                                              .naive_mul(&DensePolynomial{coeffs: vec![vec![-F::from(1)], vec![F::from(0); n-1], vec![F::from(1)]].concat()});
    let b_witness: DensePolynomial<F> =  witness.clone() + blinding_factor.clone();
    
    
        let b_witness_plus:  DensePolynomial<F> = DensePolynomial{coeffs: b_witness.coeffs.par_iter()
            .enumerate()
            .map(|(i, coeff)| coeff*g.pow(&[i as u64]))
            .collect()};
        let b_witness_plus_plus: DensePolynomial<F> = DensePolynomial{coeffs: b_witness_plus.coeffs.par_iter()
            .enumerate()
            .map(|(i, coeff)| coeff*g.pow(&[i as u64]))
            .collect()};
        
            // psi
    let psi: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("psi_coeffs.txt").unwrap() };
    //for _ in 729..1024 {
    //    coeffs.push(F::default());
    //}
    //assert_eq!(witness.evaluate(&(g*g*r)), witness_plus_plus.evaluate(&r));
    let y_start: F = b_witness.evaluate(&F::from(1));
    let y_end: F = b_witness.evaluate(&g.pow(&[728]));
    //assert_eq!(psi_challenge(&b_witness.evaluate(&r), &b_witness_plus_plus.evaluate(&r), &psi.evaluate(&r), &r, &(n as u64), &g), test.evaluate(&r));
    
    //let test: DensePolynomial<F> = mod_poly_poly(&b_witness.clone(), &b_witness_plus.clone(), n, g);
    //for i in 0..729 {
    //    if test.evaluate(&g.pow(&[i as u64])) != F::from(0) {
           // if mod_poly(b_witness.evaluate(&g.pow(&[i as u64])), b_witness_plus.evaluate(&g.pow(&[i as u64]))) != F::from(0) {
    //        println!("{}", i);
    //    }
    //}
    //return;
    //let check = test.div(&DensePolynomial{ coeffs: [vec![-F::from(1)], vec![F::from(0); n-1], vec![F::from(1)]].concat()});
    //assert_eq!(check.evaluate(&r), test.evaluate(&r)/(r.pow(&[n as u64]) - F::from(1)));
    //return;
    //let test_compare: F = mod_challenge( &b_witness.evaluate(&r), &b_witness.evaluate(&(g*r)),&r, &g, &(n as u64 ));
    //assert_eq!(test.evaluate(&r), test_compare);
    //return;
    let s_ord: u64 = 81*64;
    let rep_param: usize = 1;
    //let (paths, points, mroots, indices) = fri_prove(witness, l_list.clone(), s, r, s_ord, rep_param.clone());
    //let (b, xcc) = fri_verify(paths, points, mroots,  l_list, s, r, s_ord, rep_param as u8);
    //println!("Verification successful");
    //return;
    let now = Instant::now();
    let (challenge_vals, roots_fri, roots, paths_fri, additional_paths, points_fri, additional_points) = prove(witness, psi, g.clone(), s.clone(), r, s_ord, &y_start, &y_end, l_list.clone(), rep_param);
    println!("Prover Time: {}", now.elapsed().as_secs());
    let now = Instant::now();
    let b = verify(challenge_vals, roots_fri, roots, paths_fri, additional_paths, points_fri, additional_points, g, s, r, &729, s_ord, &y_start, &y_end, l_list, rep_param);
    println!("Verifier Time: {}", now.elapsed().as_millis());
    if b {
        println!("Verification successful");
    } else {
        println!("Verification failed");
    }
    return;
}
//println!("Hello, world!, {:?}",b);


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


fn write_to_file(interp_poly: &Vec<F>) -> std::io::Result<()> {
    let mut file = File::create("interp_poly_coeffs.txt")?;
for i in interp_poly.iter() {
    let mut i0: Fp = Fp::from(0);let mut i1: Fp = Fp::from(0);
    if !(i.c0==Fp::from(0)) {i0 = i.c0;}
    if !(i.c1==Fp::from(0)) {i1 = i.c1;}
    writeln!(file, "{}, {}", i0, i1)?;
}

Ok(())}
fn create_s() -> F {
    let mut s: F = F::new(MontFp!("2"), MontFp!("1"));

let max_pow: u64 = 18446744073709551615;
let gen_1: F = s.pow(&[max_pow]); // s^(maxx)
let gen_2: F = gen_1.pow(&[max_pow]); // s^(maxx^2)
let gen_3: F = gen_2.pow(&[max_pow]); // s^(maxx^3)
let gen_4: F = gen_3.pow(&[max_pow]); // s^(maxx^4)
let gen_5: F = gen_4.pow(&[max_pow]); // s^(maxx^5)
let gen_6: F = gen_5.pow(&[max_pow]); // s^(maxx^6)

let r_1: F = gen_1.pow(&[673548990033760218]); // s^(fac1*maxx)
let r_2: F = gen_2.pow(&[2443020786481962043]); // s^(fac2*maxx**2)
let r_3: F = gen_3.pow(&[3378466268170092361]); // s^(fac3*maxx**3)
let r_4: F = gen_4.pow(&[11301019636957581654]); // s^(fac4*maxx**4)
let r_5: F = gen_5.pow(&[7856978775279522800]); // s^(fac5*maxx**5)
let r_6: F = gen_6.pow(&[620258357900100]); // s^(fac6*maxx**6)

let mut s: F = r_1*r_2*r_3*r_4*r_5*r_6*s.pow(&[16611077425395483196]);
for j in 0..216 {
    s = s.pow(&[2]);
    if s == F::from(1){
        println!("i: {:?}", j);
        return s 
    }
}
for j in 0..131 {
    s = s.pow(&[3]);
    if s == F::from(1){
        println!("j: {:?}", j);
        return s
    }
}
s
}


use std::{
    io::{self, BufRead, BufReader},
    path::Path,
};
use std::str::FromStr;
fn lines_from_file(filename: impl AsRef<Path>) -> io::Result<Vec<F>> {
    BufReader::new(File::open(filename)?).lines()
    .map(|line| {
        let line = line?;
        let mut parts = line.trim().split(",");
        let a: Fp;
        let b: Fp;
        let a_tentative = Fp::from_str(parts.next().unwrap());
        match a_tentative {
            Ok(a_val) => a = a_val,
            Err(_) => a = Fp::from(0),
        }
        let b_tentative = Fp::from_str(parts.next().unwrap());
        match b_tentative {
            Ok(b_val) => b = b_val,
            Err(_) => b = Fp::from(0),
        }
        //let a: Fp = Fp::from_str(parts.next().unwrap()).unwrap();
        //let b: Fp = Fp::from_str(parts.next().unwrap()).unwrap();
        //println!("yes");
        Ok(F::new(a,b))
    })
    .collect()
}

fn lines_from_file_2(filename: impl AsRef<Path>) -> io::Result<Vec<F>> {
    BufReader::new(File::open(filename)?).lines()
    .map(|line| {
        let line = line?;
        let a: Fp;
        let b: Fp;
        if ! line.contains("*x") {
            a = Fp::from_str(&line).unwrap();
            b = Fp::from(0);
        }
        else if ! line.contains("+") {
            let mut parts = line.trim().split("*x");
            b = Fp::from_str(parts.next().unwrap().trim()).unwrap();
            a = Fp::from(0);
        }
        else {
        
        let mut parts = line.trim().split("*x +");
        
        b = Fp::from_str(parts.next().unwrap()).unwrap();
        //println!("b: {:?}", b);
        a = Fp::from_str(parts.next().unwrap().trim()).unwrap();
        //println!("a: {:?}", a);
        //println!("yes");
    }
        Ok(F::new(a,b))
    })
    .collect()
}
// This element has order 2^16
const FFT_GEN: F = F::new(MontFp!("17231939763216297887217622266809272467545088513556272765685947455324509609957290372982961616729012186542372231626409308935024145447"), MontFp!("14794844963276765294078403131215971792257250447333909245841623048180831260548488349233598857775988700515612163307689388891415390090"));
// This is the multiplicative generator of the field ^(l-2)
// It has order 2^217*3^136
const MULT_GEN: F = F::new(MontFp!("4887884732269044310381829002291498723817156048752319302265161467241044247866395345194043334365723689179743805338576987868462946714"), MontFp!("9775769464538088620763658004582997447634312097504638604530322934482088495732790690388086668731447378359487610677153975736925893426"));
// This function computes the FFT of a polynomial over the finite field F
use rayon::{prelude::*, join};

use crate::isogeny_prove::{initial_poly, initial_challenge};

use std::ops::Sub;
fn mod_poly_poly1(p: &DensePolynomial<F>, q: &DensePolynomial<F>, T: usize, g: F) -> DensePolynomial<F> {
    let p_squared: DensePolynomial<F> = p.naive_mul(p);
    let q_squared: DensePolynomial<F> = q.naive_mul(q);
    let p_cubed: DensePolynomial<F> = p_squared.naive_mul(p);
    let q_cubed: DensePolynomial<F> = q_squared.naive_mul(q);
    let p_squared_q: DensePolynomial<F> = p_squared.naive_mul(q);
    let q_squared_p: DensePolynomial<F> = q_squared.naive_mul(p);
    let pq: DensePolynomial<F> = p.naive_mul(q);

    let temp: DensePolynomial<F> = p_cubed+q_cubed.sub(&p_squared_q.naive_mul(q))
                                    +DensePolynomial{coeffs: vec![F::from(1488u128)]}.naive_mul(&(p_squared_q+q_squared_p))
                                    +DensePolynomial{coeffs: vec![-F::from(162000u128)]}.naive_mul(&(p_squared+q_squared))
                                    +DensePolynomial{coeffs: vec![F::from(40773375u128)]}.naive_mul(&pq)
                                    +DensePolynomial{coeffs: vec![F::from(8748000000u128)]}.naive_mul(&(p+q))
                                    +DensePolynomial{coeffs: vec![-F::from(157464000000000u128)]};

    temp.naive_mul(&DensePolynomial{coeffs: vec![-g.pow(&[T as u64-1]), F::from(1)]})
}
fn mod_poly(x: F, y: F) -> F {
    x*x*x+y*y*y-x*x*y*y+F::from(1488u128)
    *(x*x*y+y*y*x)-F::from(162000u128)*
    (x*x+y*y)+F::from(40773375u128)*x*y
    +F::from(8748000000u128)*(x+y)-
    F::from(157464000000000u128)
}