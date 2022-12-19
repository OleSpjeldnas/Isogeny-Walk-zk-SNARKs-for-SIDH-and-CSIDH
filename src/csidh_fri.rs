use ark_poly::{univariate::DensePolynomial, Polynomial};
use ark_ff::Field;
use ark_crypto_primitives::CRHScheme;
use ark_crypto_primitives::crh::poseidon;
use super::{field::F, FieldMT, poseidon_parameters};


// This function executes one folding step in the CSIDH-FRI algorithm
pub fn fold(f: &DensePolynomial<F>, l: u8, theta: F) -> DensePolynomial<F> {

let mut g_polys: Vec<Vec<F>> = Vec::new();
let d = (((Polynomial::degree(&f)+1)/(l as usize)) as f32).ceil() as usize;
for j in 0..l {
let th = theta.pow(&[j as u64]);
// The g are the g_i such that f(x) = x^i*g_i(x^l)
let mut g: Vec<F> = vec![F::from(0); d];

for (i, coeff) in f.coeffs.iter().enumerate() {
   if j == 0 && i as u8 % l == 0 {
    g[0]+= coeff*(&th);
   }
    else if (i as u8>= j) && ((i as u8 - j) % l == 0) {
        g[j as usize]+= coeff*(&th);
    }
   }
g_polys.push(g);
}
let mut final_g = vec![F::from(0); l as usize];
for j in 0..l{
for poly in g_polys.iter() {
        final_g[j as usize] += poly[j as usize];
    }
}
DensePolynomial { coeffs: final_g }
}

// Computes the Merkle tree of the folded f on the evaluation domain r<s>
pub fn round_commit(f_folded: &DensePolynomial<F>, s: &F, r: &F, s_ord: &u8) -> FieldMT {
    let leaf_crh_params = poseidon_parameters();
    let two_to_one_params = leaf_crh_params.clone();

    //let f_folded = fold(f, l, theta);
    let mut eval_vec: Vec<F> = Vec::new();
    for i in 0..*s_ord {
        eval_vec.push(f_folded.evaluate(&(r*&s.pow(&[i as u64]))));
    }
    // Let k be such that 2^k-1 < s_ord <= 2^k. Fill the 2^k-s_ord last entries of
    // eval_vec with 0 
    let k = (((*s_ord as f32).log2()).ceil()) as u8;
    for _ in 0..(2u8.pow(k.into())-s_ord) {
        eval_vec.push(F::from(0));
    }
    
    FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        eval_vec.iter().map(|x| [*x]),
    )
    .unwrap()
}

pub fn commit(f: DensePolynomial<F>, l_list: Vec<u8>, mut s: F, mut r: F, mut s_ord: u8) -> Vec<FieldMT> {
let mut MTrees: Vec<FieldMT> = Vec::new();
let params = poseidon_parameters();
let mut theta_vec: Vec<F> = Vec::new();
theta_vec.push(poseidon::CRH::<F>::evaluate(&params, f.coeffs.clone()).unwrap());
let n_rounds = l_list.len();
let mut f_folded = f;
for i in 0..n_rounds {
    f_folded = fold(&f_folded, l_list[i], theta_vec[i]);
    MTrees.push(round_commit(&f_folded, &s, &r, &s_ord));
    // raise r and s to the power of l[i]
    r = r.pow(&[l_list[i] as u64]);
    s = s.pow(&[l_list[i] as u64]);
    s_ord = s_ord/l_list[i];
    theta_vec.push(poseidon::CRH::<F>::evaluate(&params, [theta_vec.clone(),f.coeffs.clone()].concat()).unwrap());
}

MTrees
}