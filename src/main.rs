use ark_ff::{MontFp};
use ark_poly::{univariate::DensePolynomial,Polynomial};
//pub mod field;
pub mod test_field;
use ark_std::iterable::Iterable;
//use field::{FqConfig, F as Fp, Fq2 as F};
use test_field::{F as Fp, Fq2 as F};
pub mod matrix;
use matrix::*;
use merkle::{poseidon_parameters, FieldMT, FieldPath};
use rayon::prelude::*;
use std::{
    collections::hash_map::DefaultHasher,
    fs::File,
    hash::{Hash, Hasher},
    io::{self, BufRead, BufReader, Write},
    path::Path,
    str::FromStr,
    time::Instant,
};
use test_field::{Fq2 as F, F as Fp};

// TODO: Move to separate crate
pub mod csidh_fri;
pub mod get_roots;
pub mod isogeny_prove;
use isogeny_prove::{prove, verify};
use csidh_fri::*;
use ark_ff::Field;
use std::hash::{Hash, Hasher};
use ark_crypto_primitives::CRHScheme;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
use std::fs::File;
use ark_std::test_rng;
use std::io::{ Write};
//use ark_crypto_primitives::*;
pub mod merkle;
use merkle::{FieldMT, poseidon_parameters, FieldPath};
use std::time::{Instant};
use rayon::prelude::*;
use ark_serialize::{CanonicalSerialize, Compress};
use pprof::*;
fn main() {
    // l_vec contains all the folding factors
    let l_list: Vec<usize> = vec![vec![3; 6], vec![2;2]].concat();
    let mut s: F = MULT_GEN;
    
    
    for _ in 0..210 {
        s = s.pow(&[2]);
    }
    for _ in 0..131 {
        s = s.pow([3]);
    }
    let mut g: F = s.clone();
    for _ in 0..6 {
        g = g.pow(&[2]);
    }
    let r: F = F::new(Fp::from(5), Fp::from(3));

    let witness: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("new_coeffs.txt").unwrap() };
    // Witness(x+1)
    let n = witness.coeffs.len();

    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }.naive_mul(&DensePolynomial {
        coeffs: vec![vec![-F::from(1)], vec![F::from(0); n - 1], vec![F::from(1)]].concat(),
    });
    let b_witness: DensePolynomial<F> = witness.clone() + blinding_factor;

    let _b_witness_plus: DensePolynomial<F> = DensePolynomial {
        coeffs: b_witness.coeffs.par_iter().enumerate().map(|(i, coeff)| coeff * g.pow([i as u64])).collect(),
    };

    // psi
    let psi: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("psi_coeffs.txt").unwrap() };

    let y_start: F = b_witness.evaluate(&F::from(1));
    let y_end: F = b_witness.evaluate(&g.pow(&[728]));
    

    let s_ord: u64 = 729*64;
    let rep_param: usize = 2;

    //let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).build().unwrap();
    let now = Instant::now();
    let (challenge_vals, roots_fri, roots, paths_fri, additional_paths, points_fri, additional_points) =
        prove(witness, psi, g, s, r, s_ord, &y_start, &y_end, l_list.clone(), rep_param);
    println!("Prover Time: {} s", now.elapsed().as_secs());
    // if let Ok(report) = guard.unwrap().report().build() {
    //    let file = File::create("flamegraph.svg").unwrap();
    //     report.flamegraph(file).unwrap();
    //};
    let now = Instant::now();
    let b = verify(
        challenge_vals.clone(),
        roots_fri.clone(),
        roots.clone(),
        paths_fri.clone(),
        additional_paths.clone(),
        points_fri.clone(),
        additional_points.clone(),
        g,
        s,
        r,
        &729,
        s_ord,
        &y_start,
        &y_end,
        l_list,
        rep_param,
    );
    println!("Verifier Time: {} ms", now.elapsed().as_millis());
    if b {
        let size1 = challenge_vals.serialized_size(Compress::Yes);
        let size2 = roots.serialized_size(Compress::Yes);
        let size3 = points_fri.serialized_size(Compress::Yes);
        let size4 = roots_fri.serialized_size(Compress::Yes);
        let size5 = paths_fri.serialized_size(Compress::Yes);
        let size6 = additional_paths.serialized_size(Compress::Yes);
        let size7 = additional_points.serialized_size(Compress::Yes);

        println!("Proof Size: {} kB", ((size1+size2+size3+size4+size5+size6+size7) as f32)/1024f32);
        println!("Verification successful");
    } else {
        println!("Verification failed");
    }
    return;
}

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
