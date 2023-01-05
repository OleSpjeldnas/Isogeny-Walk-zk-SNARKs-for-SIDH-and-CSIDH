use ark_crypto_primitives::crh::sha256::digest::generic_array::typenum::array;
use ark_ff::{MontFp};
use ark_poly::DenseUVPolynomial;
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
use std::fs::File;
use std::io::{stdin, Write};
//use ark_crypto_primitives::*;
pub mod merkle;
use merkle::{FieldMT, poseidon_parameters, FieldPath};
use get_roots::*;

fn main() {
    let mut test_arr: Vec<(F, F)> = Vec::new();
    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();
    for i in 0..10 {
        x.push(i as f64);
        y.push(i as f64+10f64);
        test_arr.push((F::from(i), F::from(i+1)));
    }
    let lag = lagrange_interp(&x, &y);
    for z in 0..10 {

        assert!((y[z] -eval_poly(&lag, z as f64)).abs() < 1f64);
    }
    return;
    let int = lagrange_interpolation(&test_arr);
    //let int2: Vec<F> = lagrange_interpolationn(&test_arr);
    //assert_eq!(int2, int);
    println!("int: {:?}", int.len());
    let test_poly: DensePolynomial<F> = DensePolynomial { coeffs: int.clone() };
    //assert_eq!(test_poly.coeffs, int);
    for i in 0..10 {
        assert_eq!(test_poly.evaluate(&F::from(i as u64)), test_arr[i].1);
    }
    return; 
    //let interp_poly: DensePolynomial<F> = DensePolynomial{
    //    coeffs: lines_from_file("interp_poly_coeffs.txt").unwrap()
    //};
    let mut roots: Vec<F> = vec![F::from(0), F::from(54000)];
    let p: DensePolynomial<F> = DensePolynomial { coeffs: mod_poly_eval_at_y_div_by_x_minus_z(roots[0], roots[1]) };
    assert_eq!(p.evaluate(&roots[1]), F::from(0));
    let roots_arr: Vec<F> = find_roots(&mut roots, 727);
    
    
    //for i in 0..roots_arr.len()-1 {
    //    println!("i: {:?}", i);
    //    assert_eq!(mod_poly(roots_arr[i], roots_arr[i+1]), F::from(0));
    //}
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
let r: F = F::from(5);
//let mut s: F = F::new(MontFp!("3491346237335031650272735001636784802726540034823085215903686762315031605618853817995745238832659777985531289527554991334616390510"), MontFp!("10474038712005094950818205004910354408179620104469255647711060286945094816856561453987235716497979333956593868582664974003849171528"));
//let s_square: F = F::new(MontFp!("12469093704767970179545482148702802866880500124368161485370309865410827162924477921413375852973784921376897462598410683337915680391"), MontFp!("19950549927628752287272771437924484587008800198989058376592495784657323460679164674261401364758055874203035940157457093340665088625"));
//let mut s = F::from(-5).sqrt().unwrap();
//assert_eq!(F::from(7)*F::new(Fp::from(4), Fp::from(3)), F::new(Fp::from(28), Fp::from(21)));
let mut s: F = F::new(MontFp!("3"), MontFp!("1"));

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
        return;
    }
}
for j in 0..131 {
    s = s.pow(&[3]);
    if s == F::from(1){
        println!("j: {:?}", j);
        return;
    }
}
//assert_eq!(roots_arr[1], interp_poly.evaluate(&F::from(1)));
//return;
let mut interpolation_points: Vec<(F,F)> = Vec::new();
for (i, root) in roots_arr.iter().enumerate() {
        interpolation_points.push((F::from(i as u64), *root));
}
println!("here");
let mut interp_poly = lagrange_interpolation(&interpolation_points[..5]);
//interp_poly.reverse();
let pol: DensePolynomial<F> = DensePolynomial::from_coefficients_vec(interp_poly.clone());
for i in 0..5 {
    println!("i: {:?}", i);
assert_eq!(roots_arr[i], pol.evaluate(&F::from(i as u64)));}
//assert_eq!(roots_arr[i], interpolation_points[i].1);}
return;
println!("mere");
let hello = write_to_file(&interp_poly);
return;
let s_ord: u8 = 72;
//let s: F = F::new(Fp::from(17), Fp::from(0)).sqrt().unwrap();
//println!("s: {:?}", s);

//let mtrees = commit(poly, l_vec, s, r, s_ord);
//let mut mat: Array2<F> = Array::zeros((0, 10));
//let relevant_primes: Vec<u64> = vec![3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587];
//let relevant_primes: Vec<u64> = vec![4, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587];
//let tt: Vec<u64> = vec![2,5,6];
let gen = raise_to_power(s, 2u64, 213);
let gen3 = raise_to_power(gen, 3u64, 135);
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
println!("length points2:{:?}", eval2.len());
let (paths, points, roots) = prove(poly, l_vec.clone(), s, r, 72, 3);
let b = verify(paths, points, roots, l_vec, s, r, 72, 3);


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

fn lagrange_interpolation1(points: &[(F, F)]) -> Vec<F> {
    let n = points.len();
    let mut coefficients = vec![F::from(0); n];
    for i in 0..n {
        println!("i: {:?}", i);
        let mut term = F::default();
        for j in 0..n {
            if i != j {
                let mut coefficient = points[i].1 / (points[i].0 - points[j].0);
                for k in 0..n {
                    if k != i && k != j {
                        coefficient = coefficient / (points[i].0 - points[k].0);
                    }
                }
                term = term + coefficient;
            }
        }
        coefficients[i] = term;
    }
    coefficients
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

fn lagrange_interpolation2(points: &[(F, F)]) -> Vec<F> {
    let n = points.len();
    let mut w = vec![F::from(0); n];
    let mut coefficients = vec![F::from(0); n];

    for i in 0..n {
        w[i] = F::from(1) / points[i].0;
        for j in 0..i {
            w[i] = w[i] * (points[i].0 - points[j].0);
        }
        for j in (i + 1)..n {
            w[i] = w[i] * (points[i].0 - points[j].0);
        }
    }

    for i in 0..n {
        let mut numerator = F::from(0);
        let mut denominator = F::from(0);
        for j in 0..n {
            if i != j {
                let term = w[j] / (points[i].0 - points[j].0);
                numerator = numerator + term * points[j].1;
                denominator = denominator + term;
            }
        }
        coefficients[i] = numerator / denominator;
    }

    coefficients
}
use std::{
    io::{self, BufRead, BufReader},
    path::Path,
};
use std::str::FromStr;

use std::cmp::Ordering;

struct KDTree {
    point: (F, F),
    left: Option<Box<KDTree>>,
    right: Option<Box<KDTree>>,
}

impl KDTree {
    fn new(point: (F, F)) -> Self {
        KDTree {
            point,
            left: None,
            right: None,
        }
    }

    fn insert(&mut self, point: (F, F)) {
        match point.0.cmp(&self.point.0) {
            Ordering::Less => {
                if let Some(ref mut left) = self.left {
                    left.insert(point);
                } else {
                    self.left = Some(Box::new(KDTree::new(point)));
                }
            }
            Ordering::Greater => {
                if let Some(ref mut right) = self.right {
                    right.insert(point);
                } else {
                    self.right = Some(Box::new(KDTree::new(point)));
                }
            }
            Ordering::Equal => {
                self.point = point;
            }
        }
    }
}

fn lagrange_interpolation(points: &[(F, F)]) -> Vec<F> {
    let mut tree = KDTree::new(points[0]);
    for i in 1..points.len() {
        tree.insert(points[i]);
    }

    let n = points.len();
    let mut coefficients = vec![F::from(0); n];
    let mut stack = Vec::new();
    stack.push((&tree, F::from(1)));
    while let Some((node, mut term)) = stack.pop() {
        if let Some(ref left) = node.left {
            stack.push((left, term * (node.point.0 - left.point.0)));
        }
        if let Some(ref right) = node.right {
            stack.push((right, term * (node.point.0 - right.point.0)));
        }
        if node.left.is_none() && node.right.is_none() {
            let i = points.iter().position(|&p| p == node.point).unwrap();
            coefficients[i] += term * node.point.1;
        }
    }
    coefficients
}
fn lagrange_interp(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut result = vec![0.0; n];
    for i in 0..n {
        let mut c = y[i];
        for j in 0..n {
            if i == j {
                continue;
            }
            c /= (x[i] - x[j]);
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            result[j] -= c / (x[j] - x[i]);
            result[i] += c / (x[j] - x[i]);
        }
    }
    result
}
fn eval_poly(coeffs: &[f64], z: f64) -> f64 {
    let mut result = 0.0;
    let mut z_pow = 1.0;
    for &coeff in coeffs.iter().rev() {
        result += coeff * z_pow;
        z_pow *= z;
    }
    result
}
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