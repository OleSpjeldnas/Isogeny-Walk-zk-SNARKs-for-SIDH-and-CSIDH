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
use ark_crypto_primitives::CRHScheme;
use ark_sponge::poseidon::PoseidonConfig;
use ark_crypto_primitives::crh::poseidon;
use ark_ff::UniformRand;
use std::fs::File;
use std::io::{stdin, Write};
//use ark_crypto_primitives::*;
pub mod merkle;
use merkle::{FieldMT, poseidon_parameters, FieldPath};
use get_roots::*;
use std::time::{Duration, Instant};

fn main() {
    // l_vec contains all the folding factors
    let l_list: Vec<usize> = vec![vec![3; 4], vec![4; 2], vec![2]].concat();
    let mut s: F = MULT_GEN;
    
    for i in 0..130 {
        s = s.pow(&[3]);
    }
    let g: F = s.clone();
    for i in 0..212 {
        s = s.pow(&[2]);
    }
    let r: F = F::new(Fp::from(5), Fp::from(3));
    //let order_group: u64 = 2u64.pow(6);
    //let w = FFT_GEN.pow(&[order_group]).inverse().unwrap();
    // Witness polynomial
    let witness: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("new_coeffs.txt").unwrap() };
    // Witness(x+1)
    let witness_plus: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("coeffs_plus.txt").unwrap() };
    // Witness(x+2)
    let witness_plus_plus: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("coeffs_plus_plus.txt").unwrap() };
    // psi
    let psi: DensePolynomial<F> = DensePolynomial { coeffs: lines_from_file_2("psi_coeffs.txt").unwrap() };
    //for _ in 729..1024 {
    //    coeffs.push(F::default());
    //}
    let y_start: F = witness.evaluate(&F::from(1));
    let y_end: F = witness.evaluate(&g.pow(&[728]));
    let s_ord: u64 = 729*32;
    let rep_param: usize = 3;
    let now = Instant::now();
    let (challenge_vals, roots_fri, roots, paths_fri, additional_paths, points_fri, additional_points) = prove(witness, witness_plus, witness_plus_plus, psi, g, s, r, s_ord, &y_start, &y_end, l_list.clone(), rep_param);
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
    //let mut eval_vec: Vec<F> = Vec::new();
    //let w = FFT_GEN.pow(&[2u64.pow(6)]);
    //for i in 0..1024 {
    //    eval_vec.push(witness.evaluate(&w.pow(&[i as u64])));
    //}
    //let check = fft(&coeffs);
    //let ver = ifft(&check);
    //let mut ver2: Vec<F> = Vec::new();
    //let n: u64 = check.len() as u64;
    //let n_inv: F = F::from(1024).inverse().unwrap();
    //for coeff in ver.iter() {
    //    ver2.push(*coeff* n_inv);
    //}
    
    //assert_eq!(coeffs, ver2);
    return;
    // Polynomial p s.t. p(s^i) = witness(s^(i+1))
    let coeffs_plus: Vec<F> = lines_from_file_2("coeffs_plus.txt").unwrap();
    let witness_plus: DensePolynomial<F> = DensePolynomial { coeffs: coeffs_plus.clone() };
    //let test = fft_multiply(&coeffs, &coeffs_plus);
    //let test_poly = DensePolynomial { coeffs: test.clone() };

    let test_point: F = F::from(243);
    //let test_1: Vec<F> = vec![F::from(2), F::from(3), F::from(1)];
    //let test_2: Vec<F> = vec![F::from(1), F::from(2), F::from(1)];
    //let poly_mult = tom_cook(&coeffs, &coeffs_plus);
    //let poly_mult2 = toom_cook_mul(&coeffs, &coeffs_plus);
    //let test_poly = DensePolynomial { coeffs: poly_mult2.clone() };
    //let poly_test = witness.naive_mul(&witness_plus);
    //let t = naive_mul2(&coeffs, &coeffs_plus);
    //assert_eq!(poly_mult2, t);
    //assert_eq!(witness.naive_mul(&witness_plus).coeffs[0], test[0]);
    //assert_eq!(test_poly.evaluate(&test_point), witness.evaluate(&test_point)*(&witness_plus.evaluate(&test_point)));
    return;
    //let karatsuba_test = karatsuba(&coeffs, &coeffs_plus);
    //let c_2 = mod_poly_quotient(witness.clone(), witness_plus.clone());
    return;
    let gen: F = F::new(MontFp!("12506453529958197004077348547101724167733006220309459995744807787144199909522473787798149509471425651348760555246425049364872841517"), MontFp!("21367825699298292014944120575290706677154228414033088926289743967671269667328826095717704544576378050739679286210174567998522892433"));
    for i in 0..728 {
        println!("i:{}", i);
        assert_eq!(mod_poly(witness.evaluate(&gen.pow(&[i as u64])),witness_plus.evaluate(&gen.pow(&[i as u64]))), F::from(0));
    }
    return;

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
let l_vec: Vec<usize> = vec![2,3,2];
let r: F = F::from(3);
let s_ord: u8 = 72;

//let mtrees = commit(poly, l_vec, s, r, s_ord);
//let mut mat: Array2<F> = Array::zeros((0, 10));
//let relevant_primes: Vec<u64> = vec![3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587];
//let relevant_primes: Vec<u64> = vec![4, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587];
//let tt: Vec<u64> = vec![2,5,6];
//let gen = raise_to_power(s, 2u64, 213);
//let gen3 = raise_to_power(gen, 3u64, 135);
//let (mt, f1, )
let f1 = fold(&poly, 2, F::from(2));
//let (mt1, r1, eval1) = round_commit(&poly, &gen3, &r, &s_ord);
//let (mt2, r2, eval2) = round_commit(&f1, &gen3.pow(&[2]), &r.pow(&[2]), &(s_ord/2));
let index: usize = 5;
//println!("length points1:{:?}", eval1.len());
//println!("length points2:{:?}", eval2.len());
//let mut points: Vec<F> = Vec::new();
//points.push(eval2[5]);
//points.push(eval1[5]);
//points.push(eval1[17]);
//let t: F = gen3.pow(&[12]);
//assert_eq!(points[0], f1.evaluate(&(r.pow(&[2])*gen3.pow(&[10]))));
//assert_eq!(points[1], poly.evaluate(&(r*gen3.pow(&[5]))));
//assert_eq!(points[2], poly.evaluate(&(r*gen3.pow(&[5])*t)));
//println!("length points2:{:?}", eval2.len());
//let (paths, points, roots) = prove(poly, l_vec.clone(), s, r, 72, 3);
//let b = verify(paths, points, roots, l_vec, s, r, 72, 3);


//println!("Hello, world!, {:?}",b);

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

fn mod_poly(x: F, y: F) -> F {
    x*x*x+y*y*y-x*x*y*y+F::from(1488u128)
    *(x*x*y+y*y*x)-F::from(162000u128)*
    (x*x+y*y)+F::from(40773375u128)*x*y
    +F::from(8748000000u128)*(x+y)-
    F::from(157464000000000u128)
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
fn fft(a: &[F]) -> Vec<F> {
    let n = a.len();
    //println!("h_hat: {}", n);
    if n == 1 {
        return a.to_vec();
    }
    let mut a0 = vec![F::from(0); n/2];
    let mut a1 = vec![F::from(0); n/2];
    for i in 0..n/2 {
        a0[i] = a[i*2];
        a1[i] = a[i*2+1];
    }
    let y0 = fft(&a0);
    let y1 = fft(&a1);
    let mut y = vec![F::from(0); n];
    let order_group: u64 = 2u64.pow(16) / (n as u64);
    let mut w: F = FFT_GEN.pow(&[order_group]);
    let mut wk = F::from(1);
    for k in 0..n/2 {
        y[k] = y0[k] + wk * y1[k];
        y[k+n/2] = y0[k] - wk * y1[k];
        wk = wk * w;
    }
    y
}

fn ifft(a: &[F]) -> Vec<F> {
    let n = a.len();
    if n == 1 {
        return a.to_vec();
    }
    let mut a0 = vec![F::from(0); n/2];
    let mut a1 = vec![F::from(0); n/2];
    for i in 0..n/2 {
        a0[i] = a[i*2];
        a1[i] = a[i*2+1];
    }
    let y0 = ifft(&a0);
    let y1 = ifft(&a1);
    let mut y = vec![F::from(0); n];
    let order_group: u64 = 2u64.pow(16) / (n as u64);
    let w = FFT_GEN.pow(&[order_group]).inverse().unwrap();
    let mut wk = F::from(1);
    for k in 0..n/2 {
        y[k] = y0[k] + wk * y1[k];
        y[k+n/2] = y0[k] - wk * y1[k];
        wk = wk * w;
    }
    y
}

// Multiply two polynomials f and g
fn fft_multiply(f: &[F], g: &[F]) -> Vec<F> {
    let m = f.len();
    let n = g.len();
    let mut f = f.to_vec();
    let mut g = g.to_vec();
    // Pad the input polynomials with zeros to the nearest power of 2

    
    let k1 = (((m as f32).log2()).ceil()) as u8 + 1;
    println!("len_before: {}", f.len());
    for _ in 0..2usize.pow(k1 as u32) - m {
        f.push(F::from(0));
    }
    for _ in 0..2usize.pow(k1 as u32) - n {
        g.push(F::from(0));
    }
    //f.resize(2usize.pow(k1 as u32), F::from(0));
    //g.resize(2usize.pow(k2 as u32), F::from(0));
    // Compute the FFTs of f and g
    let f_hat = fft(&f);
    let g_hat = fft(&g);
    // Multiply the FFTs pointwise
    
    let h_hat: Vec<F> = f_hat.iter().zip(g_hat.iter()).map(|(f, g)| *f * *g).collect();
    // Compute the inverse FFT of the result
    println!("len_after: {}", h_hat.len());
    
    let h = ifft(&h_hat);
    // Truncate the result to the correct degree
    h[..m+n-1].to_vec()
}
