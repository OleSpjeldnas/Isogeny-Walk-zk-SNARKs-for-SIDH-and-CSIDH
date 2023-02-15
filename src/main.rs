use ark_ff::MontFp;
use ark_poly::univariate::DensePolynomial;
//pub mod field;
pub mod sidh_sike_p434;
use sidh_sike_p434::{F as Fp, Fq2 as F};
pub mod matrix;
use matrix::*;
pub mod run_isogeny;
use run_isogeny::*;
use merkle::{poseidon_parameters, FieldMT, FieldPath};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, Result},
    path::Path,
    str::FromStr,
    time::Instant,
};
use std::ops::Div;
pub mod csidh_512;
use csidh_512::{K, Kp};
// TODO: Move to separate crate
pub mod generalized_fri;
pub mod get_roots;
//pub mod isogeny_prove;
pub mod csidh_prove;
use csidh_prove::{prove, verify};
//use isogeny_prove::{prove, verify};
use generalized_fri::*;
use ark_ff::Field;
use ark_crypto_primitives::CRHScheme;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
use ark_poly::Polynomial;
use ark_std::test_rng;
pub mod merkle;
use ark_serialize::{CanonicalSerialize, Compress};

fn load_mod_polys(filename: impl AsRef<Path>) 
-> Result<Vec<(u32, Vec<(u32, K)>)>> 
{
    BufReader::new(File::open(filename)?).lines()
    .map(|line| {
        //let line = line?.replace("'", "").replace("[", "");
        let line = line?.replace("'", "");
        let line = line.split_once(",").unwrap();
        let x_pow: u32 = u32::from_str(line.0.replace("[", "").trim()).unwrap();
        let y_coeff: Vec<(u32, K)> = line.1.split("], [").map(|x| {
            let x = x.replace("]", "").replace("[", "");
            let mut parts = x.split(", ");
            let y: u32 = u32::from_str(parts.next().unwrap().trim()).unwrap();
            let pre_val = parts.next().unwrap().trim();
            let val: Kp;
            if pre_val.chars().next().unwrap() == char::from_str("-").unwrap(){
                val = -Kp::from_str(&pre_val[1..]).unwrap();
                
            }
            else {
                val = Kp::from_str(pre_val).unwrap();}
                //println!("val: {:?}", val);
            (y, K::new(val, Kp::from(0)))
        }).collect();
        Ok((x_pow, y_coeff))
    })
    .collect()
}
fn main() {
    let mut witnesses: Vec<DensePolynomial<K>> = Vec::new();
    let mut witnesses_g: Vec<DensePolynomial<K>> = Vec::new();
    for i in 0..589 {
        let filename = format!("../Witness_Powers/{}.txt", i);
        witnesses.push(DensePolynomial{coeffs: load_witness(filename).unwrap()});
        let filename = format!("../Witness_Powers_g/{}.txt", i);
        witnesses_g.push(DensePolynomial{coeffs: load_witness(filename).unwrap()});
    }
    let primes: Vec<u32> = vec![3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,587];
    
    let mod_polys: Vec<Vec<(u32, Vec<(u32, K)>)>> = primes.iter().map(|x| {
        let filename = format!("../Mod_Polys_Univariate/{}.txt", x);
        load_mod_polys(filename).unwrap()
    }).collect();
let g: K = K::new(MontFp!("359043299558979565001571415256768305303404570834255644915881154600034969116473755419281641131027150698946785484257017053947786613776504114272667938230180"), MontFp!("4823725014927225578723092749633653720120879328630769151961249797747321590005955903193275095621673255443071378032591992903550006245259265723755451585798774"));
// Point of order 510510*19 = 9699690/ 74*5*3
let s: K = K::new(MontFp!("1836047270781180842167226749359467367114963289807460284073448736347315926240287364897381624845231345995631182615893220800397820653140409209944876401904601"), MontFp!("1557093291095712624420063998101784186513188938002538018122484637471338545927556000331103701274394862447229633224212846627286518669976988296359081806915762"));
let r: K = K::from(2);
let w_r:K = witnesses[1].evaluate(&r);
let w_gr:K = witnesses[1].evaluate(&(g*r));
let quotient_poly: DensePolynomial<K> = transition_constraint(&witnesses, &witnesses_g, &g, &5);
let v_schek: K = quotient_poly.evaluate(&r);
let vv: K = transition_verifier(&w_r, &w_gr, &r, &g, &5);
assert_eq!(v_schek, vv);
return;
let t: u8 = 5;
//let s_ord: u64 = 9699690;
let s_ord: u64 = 74*5*3;
let y_start: K = witnesses[1].evaluate(&K::from(1));
let y_end: K = witnesses[1].evaluate(&K::from(369));
//let l_list: Vec<usize> = vec![2, 3, 5, 7, 11, 13, 17];
let l_list: Vec<usize> = vec![3, 5];
let rep_param: usize = 1;
let grinding_param: u8 = 3;

let now = Instant::now();
    let (challenge_vals, 
        roots_fri, 
        roots, 
        paths_fri, 
        additional_paths, 
        points_fri, 
        additional_points) =
        prove(&mod_polys, witnesses, witnesses_g, g, &t, s, r, &s_ord, &y_start, &y_end, &l_list, rep_param, grinding_param);
    println!("Prover Time: {} s", now.elapsed().as_secs());
    let now = Instant::now();
    let b = verify(&mod_polys,
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
        s_ord,
        5,
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

}

fn load_witness(filename: impl AsRef<Path>) -> io::Result<Vec<K>> {
    BufReader::new(File::open(filename)?).lines()
    .map(|line| {
        let line = line?;
        let mut parts = line.trim().split(",");
        let a_str = parts.next().unwrap().trim();
        let b_str =  parts.next().unwrap().trim();
        let a = Kp::from_str(a_str).unwrap();
        let b = Kp::from_str(b_str).unwrap();
        Ok(K::new(a,b))
    }).collect()
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

fn transition_constraint(witnesses: &Vec<DensePolynomial<K>>, witnesses_g: &Vec<DensePolynomial<K>>, g: &K, t: &usize)
 -> DensePolynomial<K> 
{   let mut g_pow = K::from(1);
    //
    let primes: Vec<u32> = vec![3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,587];
    //let primes: Vec<u32> = vec![3,5];
    let mod_polys: Vec<Vec<(u32, Vec<(u32, K)>)>> = primes.iter().map(|x| {
        let filename = format!("../Mod_Polys_Univariate/{}.txt", x);
        load_mod_polys(filename).unwrap()
    }).collect();
    let mut final_poly = DensePolynomial{coeffs: vec![K::from(0)]};
    for mod_poly in mod_polys.iter() {
        let mut temp_poly: DensePolynomial<K> = DensePolynomial{
            coeffs: vec![K::from(0)]
        };
        for (x_pow, y_term) in mod_poly.iter() {
            let mut f_y: DensePolynomial<K> = DensePolynomial{
                coeffs: vec![K::from(0)]
            };
            for (y_pow, coeff) in y_term.iter() {
                f_y = f_y + witnesses_g.clone()[*y_pow as usize].naive_mul(&DensePolynomial{
                    coeffs: vec![*coeff]
                });
            }
            temp_poly = temp_poly + witnesses.clone()[*x_pow as usize].naive_mul(&f_y);
            
        }
        let mut div_poly = DensePolynomial{
            coeffs: vec![K::from(1)]
        };
        if final_poly.coeffs[0] == K::from(0){
            for _ in 0..*t-1 {
                div_poly = div_poly.naive_mul(&DensePolynomial{
                    coeffs: vec![-g_pow, K::from(1)]
    
                });
                g_pow = g_pow * g;
            }
        }
        else {
            for _ in 0..*t {
                div_poly = div_poly.naive_mul(&DensePolynomial{
                    coeffs: vec![-g_pow, K::from(1)]
    
                });
                g_pow = g_pow * g;
            }
        }
        final_poly = final_poly + temp_poly.div(&div_poly);
    }
    final_poly
}

fn transition_verifier(y: &K, y_g: &K, z: &K, g: &K, t: &usize)
 -> K 
{   let mut ys: Vec<K> = vec![K::from(1)];
    let mut y_gs: Vec<K> = vec![K::from(1)];
    for _ in 0..588 {
        ys.push(ys.last().unwrap() * y);
        y_gs.push(y_gs.last().unwrap() * y_g);
    }
    let mut g_pow = K::from(1);
    let primes: Vec<u32> = vec![3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,587];
    //let primes: Vec<u32> = vec![3,5];
    let mod_polys: Vec<Vec<(u32, Vec<(u32, K)>)>> = primes.iter().map(|x| {
        let filename = format!("../Mod_Polys_Univariate/{}.txt", x);
        load_mod_polys(filename).unwrap()
    }).collect();
    let mut final_val: K = K::from(0);
    for mod_poly in mod_polys.iter() {
        let mut temp_val: K = K::from(0);
        for (x_pow, y_term) in mod_poly.iter() {
            let mut f_y: K = K::from(0);
            for (y_pow, coeff) in y_term.iter() {
                f_y = f_y + y_gs[*y_pow as usize] * coeff;
            }
            temp_val = temp_val + ys[*x_pow as usize] * f_y;
            
        }
        let mut div_val = K::from(1);
        if final_val == K::from(0) {
            for _ in 0..*t-1 {
                div_val = div_val * (z-g_pow);
                g_pow = g_pow * g;
            }
        }
        else {
            for _ in 0..*t{
                div_val = div_val * (z-g_pow);
                g_pow = g_pow * g;
            }
        }
        final_val = final_val + temp_val / div_val;
    }
    final_val
}

fn eval_mod_poly(y: &K, y_g: &K, l: &u32)
 -> K 
{   let mut ys: Vec<K> = vec![K::from(1)];
    let mut y_gs: Vec<K> = vec![K::from(1)];
    for _ in 0..l+1 {
        ys.push(ys.last().unwrap() * y);
        y_gs.push(y_gs.last().unwrap() * y_g);
    }
    //let primes: Vec<u32> = vec![2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,587];
   let filename = format!("Mod_Polys_Univariate/{}.txt", l);
   let mod_poly = load_mod_polys(filename).unwrap();
    //let mut final_val: K = K::from(0);
        let mut temp_val: K = K::from(0);
        for (x_pow, y_term) in mod_poly.iter() {
            let mut f_y: K = K::from(0);
            for (y_pow, coeff) in y_term.iter() {
                f_y = f_y + y_gs[*y_pow as usize] * coeff;
            }
            temp_val = temp_val + ys[*x_pow as usize] * f_y;
            
        }
       temp_val
}