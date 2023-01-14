use ark_ff::MontFp;
use ark_poly::univariate::DensePolynomial;
//pub mod field;
pub mod sidh_sike_p434;
use sidh_sike_p434::{F as Fp, Fq2 as F};
pub mod matrix;
use matrix::*;
use merkle::{poseidon_parameters, FieldMT, FieldPath};
use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
    str::FromStr,
    time::Instant,
};

// TODO: Move to separate crate
pub mod generalized_fri;
pub mod get_roots;
pub mod isogeny_prove;
use isogeny_prove::{prove, verify};
use generalized_fri::*;
use ark_ff::Field;
use ark_crypto_primitives::CRHScheme;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
use ark_poly::Polynomial;
use ark_std::test_rng;
pub mod merkle;
use ark_serialize::{CanonicalSerialize, Compress};
fn main() {
    // l_vec contains all the folding factors
    let l_list: Vec<usize> = vec![vec![3; 6], vec![2;2]].concat();
    let mut s: F = MULT_GEN;
    
    
    for _ in 0..211 {
        s = s.pow(&[2]);
    }
    for _ in 0..131 {
        s = s.pow([3]);
    }
    let mut g: F = s.clone();
    for _ in 0..5 {
        g = g.pow(&[2]);
    }
    //for _ in 0..2 {
    //    g = g.pow(&[3]);
    //}
    let r: F = F::new(Fp::from(5), Fp::from(3));

    let witness: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file("new_coeffs.txt").unwrap() };
    let n = witness.coeffs.len();

    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }.naive_mul(&DensePolynomial {
        coeffs: vec![vec![-F::from(1)], vec![F::from(0); n - 1], vec![F::from(1)]].concat(),
    });
    let b_witness: DensePolynomial<F> = witness.clone() + blinding_factor;

    // psi
    let psi: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file("psi_coeffs.txt").unwrap() };

    let y_start: F = b_witness.evaluate(&F::from(1));
    let y_end: F = b_witness.evaluate(&g.pow(&[728]));
    

    let s_ord: u64 = 729*32;
    let rep_param: usize = 1;
    let grinding_param: u8 = 32;

    let now = Instant::now();
    let (challenge_vals, 
        roots_fri, 
        roots, 
        paths_fri, 
        additional_paths, 
        points_fri, 
        additional_points) =
        prove(witness, psi, g, s, r, s_ord, &y_start, &y_end, l_list.clone(), rep_param, grinding_param);
    println!("Prover Time: {} s", now.elapsed().as_secs());
    
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
        grinding_param
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

fn lines_from_file(filename: impl AsRef<Path>) -> io::Result<Vec<F>> {
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
// This is the multiplicative generator of the field ^(l-2)
// It has order 2^217*3^136
const MULT_GEN: F = F::new(MontFp!("4887884732269044310381829002291498723817156048752319302265161467241044247866395345194043334365723689179743805338576987868462946714"), MontFp!("9775769464538088620763658004582997447634312097504638604530322934482088495732790690388086668731447378359487610677153975736925893426"));

