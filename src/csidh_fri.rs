use ark_poly::{univariate::DensePolynomial, Polynomial};
use ark_ff::Field;
use ark_crypto_primitives::{CRHScheme, MerkleTree};
use ark_crypto_primitives::crh::poseidon;
use super::{F, FieldMT, poseidon_parameters, FieldPath, solve_linear_system, Fp};
use ndarray::prelude::*;
use ndarray::{Array, ArrayView, Axis};
use ndarray_linalg::*;


// This function executes one folding step in the CSIDH-FRI algorithm
pub fn fold(f: &DensePolynomial<F>, l: u8, theta: F) -> DensePolynomial<F> {

let mut g_polys: Vec<Vec<F>> = Vec::new();
let d = (((Polynomial::degree(f))/(l as usize)) as f32).floor() as usize + 1;
//println!("deg: {:?}", Polynomial::degree(f));
//println!("d: {:?}", d);
for j in 0..l {
let th = theta.pow(&[j as u64]);
// The g are the g_i such that f(x) = x^i*g_i(x^l)
let mut g: Vec<F> = vec![F::from(0); d];

for (i, coeff) in f.coeffs.iter().enumerate() {
    //println!("Coeff: {}", coeff);
   if (j > i as u8) && ((j - i as u8) % l == 0){
        g[i]+= coeff*(&th);
   }
    else if (i as u8>= j) && ((i as u8 - j) % l == 0) {
        g[(i-j as usize)/(l as usize)]+= coeff*(&th);
    }
   }
   //println!("Length g: {:?}", g.len());
g_polys.push(g);
}
let mut final_g = vec![F::from(0); d as usize];
for j in 0..d{
for poly in g_polys.iter() {
        final_g[j as usize] += poly[j as usize];
    }
}
while final_g.last().unwrap() == &F::from(0) {
    final_g.remove(final_g.len()-1);
}
//println!("Length final_g: {:?}", final_g.len());
DensePolynomial { coeffs: final_g }
}

// Computes the Merkle tree of the folded f on the evaluation domain r<s>
pub fn round_commit(f_folded: &DensePolynomial<F>, s: &F, r: &F, s_ord: &u8) -> (FieldMT, Fp, Vec<F>) {
    let leaf_crh_params = poseidon_parameters();
    let two_to_one_params = leaf_crh_params.clone();

    //let f_folded = fold(f, l, theta);
    let mut eval_vec: Vec<Vec<Fp>> = Vec::new();
    let mut points_vec: Vec<F> = Vec::new();
    for i in 0..*s_ord {
        let temp: F = f_folded.evaluate(&(r*&s.pow(&[i as u64])));
        eval_vec.push(vec![temp.c0, temp.c1]);
        points_vec.push(temp);
    }
    // Let k be such that 2^k-1 < s_ord <= 2^k. Fill the 2^k-s_ord last entries of
    // eval_vec with 0 
    let k = (((*s_ord as f32).log2()).ceil()) as u8;
    //println!("k:{:?}", k);
    //println!("s_ord:{:?}", s_ord);
    for _ in 0..(2u8.pow(k.into())-s_ord) {
        eval_vec.push(vec![Fp::from(0),Fp::from(0)]);
    }
    let mtree: FieldMT = FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        eval_vec.iter().map(|x| x.as_slice()),
    )
    .unwrap();
    (mtree.clone(), mtree.root(), points_vec)
}

// Returns (Merkle Trees, Merkle Roots, Evaluations)
pub fn commit(f: DensePolynomial<F>, l_list: Vec<usize>, mut s: F, mut r: F, mut s_ord: u8) -> (Vec<FieldMT>, Vec<Fp>, Vec<Vec<F>>) {
let mut mtrees: Vec<FieldMT> = Vec::new();
let mut points: Vec<Vec<F>> = Vec::new();
let mut roots: Vec<Fp> = Vec::new();

let (first_mt, first_root, first_points) = round_commit(&f, &s, &r, &s_ord);
mtrees.push(first_mt);
points.push(first_points);
roots.push(first_root);
let params = poseidon_parameters();
let mut theta_vec: Vec<F> = Vec::new();
theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0)));
let n_rounds = l_list.len();
let folded_polys: &mut Vec<DensePolynomial<F>> = &mut Vec::new();
folded_polys.push(f);
for i in 0..n_rounds {
    folded_polys.push(fold(folded_polys.last().unwrap(), l_list[i] as u8, theta_vec[i]));
    r = r.pow(&[l_list[i] as u64]);
    s = s.pow(&[l_list[i] as u64]);
    s_ord = s_ord/(l_list[i] as u8);
   let (m, r, p) = round_commit(folded_polys.last().unwrap(), &s, &r, &s_ord);
    roots.push(r);
    mtrees.push(m);
    points.push(p);
    theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0)));
}
(mtrees, roots, points)
}

pub fn query_at_index(mt1: FieldMT, mt2: FieldMT, points1: Vec<F>, points2: Vec<F>, index: usize,l: usize, n: usize) -> (Vec<FieldPath>, Vec<F>) {
    let mut paths: Vec<FieldPath> = Vec::new();
    // Query MT2 at position index, save point to var, path to vec
    let path_main: FieldPath = mt2.generate_proof(index % (n/l)).unwrap();
    paths.push(path_main);
    let y: F = points2[index % (n/l)];
    // Query MT1 at positions (index +kn/l) % n, save points to vec, paths to vec paths
    let mut points: Vec<F> = vec![y];
    for k in 0..l {
        let i: usize = (index + (k*n)/l) % n;
        let path_current: FieldPath = mt1.generate_proof(i).unwrap();
        paths.push(path_current);
        points.push(points1[i]);
    }

    // Output paths to be verified later
    (paths, points)
}

pub fn query(mtrees: Vec<FieldMT>, points: Vec<Vec<F>>, l_list: Vec<usize>, n: usize, alpha: usize) -> (Vec<FieldPath>, Vec<F>) {
let mut paths: Vec<FieldPath> = Vec::new();
let mut queried_points: Vec<F> = Vec::new();
//let n = l_list.len();
let mut indices: Vec<u64> = Vec::new();
indices.push(calculate_hash(&l_list, n as u64));
for _ in 0..alpha {
    let mut s_ord = n.clone();
    let index = *indices.last().unwrap() as usize;
    let (m, p) = query_at_index(mtrees[0].clone(), mtrees[1].clone(), points[0].clone(), points[1].clone(), index, l_list[0], s_ord);
    paths = [paths.clone(),m].concat();
    queried_points = [queried_points.clone(),p].concat();
    indices.push(calculate_hash(&indices, s_ord as u64));
    s_ord /= l_list[0];
    for (i, l) in l_list[1..].to_vec().iter().enumerate() {
        //println!("s_ord:{:?}", s_ord);
        let index = *indices.last().unwrap() as usize;
        let (m, p) = query_at_index(mtrees[i].clone(), mtrees[i+1].clone(), points[i].clone(), points[i+1].clone(), index, l_list[i], s_ord);
        paths = [paths.clone(),m].concat();
        queried_points = [queried_points.clone(),p].concat();
        indices.push(calculate_hash(&indices, s_ord as u64));
        s_ord /= l;
    }
}
(paths, queried_points)
}

pub fn prove(f: DensePolynomial<F>, l_list: Vec<usize>, s: F, r: F, s_ord: u8, alpha: usize) -> (Vec<FieldPath>, Vec<F>, Vec<Fp>) {
    let (mtrees, mroots, evals) = commit(f, l_list.clone(), s, r, s_ord.clone());

    let (paths, points) = query(mtrees, evals, l_list, s_ord as usize, alpha);

    (paths, points, mroots)
}  

    // Assert the equality y != \sum_i t^i*x^i*v_i
pub fn verify_fold_at_index(points: Vec<F>, x: F, t: F, l: usize, theta: F) -> bool {
    let y: F = points[0];
    let z_vec: Vec<F> = points[1..].to_vec();
    assert_eq!(z_vec.len(), l);
    println!("here");
    let mut mat: Vec<Vec<F>> = Vec::new();
    for i in 0..l {
        let mut vec_i: Vec<F> = Vec::new();
        let z: F = t.pow(&[i as u64])*x;
    for j in 0..l {
        vec_i.push(z.pow(&[j as u64]));
        }
        mat.push(vec_i);
}
let g_vec = solve_linear_system(mat, z_vec);
let mut y_supposedly: F = F::from(0);
for (i, val) in g_vec.iter().enumerate() {
    y_supposedly += theta.pow(&[i as u64])*val;
}
println!("Res: {:?}", y == y_supposedly);
y == y_supposedly
}

pub fn verify(mut paths: Vec<FieldPath>, mut queried_points: Vec<F>, roots: Vec<Fp>, l_list: Vec<usize>, s: F, r: F, s_ord: usize, alpha: u8) -> bool {
    let leaf_crh_params = poseidon_parameters();
    let two_to_one_params = leaf_crh_params.clone();
    //let n = l_list.len();

    let mut t_vals: Vec<F> = Vec::new();
    let mut s_vals: Vec<F> = Vec::new();
    let mut r_vals: Vec<F> = Vec::new();
    let mut s_ord_vals: Vec<usize> = Vec::new();

    s_vals.push(s); s_ord_vals.push(s_ord); r_vals.push(r);
    for (i, l) in l_list.as_slice().iter().enumerate(){
        let t_ord = s_ord/l;
        t_vals.push(s.pow(&[t_ord as u64]));
        //assert_eq!(s.pow(&[t_ord as u64]), F::from(1));
        //assert_eq!(s, F::from(1));
        
        r_vals.push(r_vals[i].pow(&[*l as u64]));
        s_vals.push(s_vals[i].pow(&[*l as u64]));
        s_ord_vals.push(s_ord_vals[i]/l);
    }
    // Define all the thetas using Fiat-Shamir
    let mut theta_vec: Vec<F> = Vec::new();
    let params = poseidon_parameters();
    let mut rr: Vec<Fp> = vec![roots[0]];
    theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, rr[..1].to_vec()).unwrap(), Fp::from(0)));
    for root in roots[1..].to_vec().iter() {
        rr.push(*root);
        theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, rr.clone()).unwrap(), Fp::from(0)));
    }
    let mut indices: Vec<u64> = vec![calculate_hash(&l_list, s_ord as u64)];
    let mut i: usize = 0;
    assert_eq!(paths.len(), queried_points.len());
    for __ in 0..alpha { 
        for (j, l) in l_list.iter().enumerate() { 

            println!("l: {:?}", l);
            println!("j: {:?}", j);
            assert!(
                paths[i].verify(
                    &leaf_crh_params,
                    &two_to_one_params,
                    &roots[j+1],
                    [queried_points[i].c0, queried_points[i].c1]
                ).unwrap());
                i+=1;
                println!("one");
            for _ in 0..*l {
    assert!(
        paths[i].verify(
            &leaf_crh_params,
            &two_to_one_params,
            &roots[j],
            [queried_points[i].c0, queried_points[i].c1]
        ).unwrap());
        i+=1;
    }
}
}

    for _ in 0..alpha {
        for (i, l) in l_list.as_slice().iter().enumerate() {
            let index = *indices.last().unwrap();
            assert!(verify_fold_at_index(
            queried_points[0..*l as usize+1].to_vec(),
            s_vals[i].pow(&[index])*r_vals[i],
            t_vals[i],
            *l as usize,
            theta_vec[i]
        ) 
    );
        indices.push(calculate_hash(&indices, s_ord as u64));
    for j in 1..l+1{
        assert!(
            paths[j as usize].verify(
                &leaf_crh_params,
                &two_to_one_params,
                &roots[i+1],
                [queried_points[j as usize].c0, queried_points[j as usize].c1]
            ).unwrap());
    //qp.drain(0..l+1);
        }
        paths.drain(0..*l as usize+1);
        queried_points.drain(0..*l as usize+1);
    }
    
    }
    return true
}

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
fn calculate_hash<T: Hash>(t: &T, n: u64) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish() % n
}