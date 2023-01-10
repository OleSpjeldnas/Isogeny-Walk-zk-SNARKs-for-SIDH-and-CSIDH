use std::ops::{Div, Mul, Sub};
use rayon::prelude::*;

use super::*;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::{UniformRand, FftField};
use ark_std::test_rng;
use merkle::{FieldMT, poseidon_parameters};
// Witness is the witness polynomial, psi the inverse of w(x)-w(g^2x), g the generator of the interpolation domain, 
//the evaluation domain is E = r<s>. Finally, s_ord is the size of E.   ->    (Challenges, roots, roots_fri, paths_fri, paths, points_fri, points)
pub fn prove(witness: DensePolynomial<F>, psi: DensePolynomial<F>, g: F, s: F, r: F, s_ord: u64, y_start: &F, y_end: &F, l_list: Vec<usize>, rep_param: usize)
 -> (Vec<F>, Vec<Fp>, Vec<Fp>, Vec<FieldPath>, Vec<Vec<FieldPath>>, Vec<F>, Vec<Vec<F>>){
   
    let n: usize = witness.coeffs.len();
   
    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    //Blind the witness
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
        // Commit to the blinded witness and psi
    let params = poseidon_parameters();
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();


    let D_0: Vec<F> = (0..s_ord).into_par_iter()
                                .map(|i| r*s.pow([i]))
                                .collect();
    let witness_evals: Vec<F> = D_0.clone().into_par_iter()
                                   .map(|x| b_witness.evaluate(&x))
                                   .collect();
    let psi_evals: Vec<F> = D_0.clone().into_par_iter()
    .map(|x| psi.evaluate(&x))
    .collect();
    println!("Checkpoint 0");

    let mut witness_evals_merkle: Vec<Vec<Fp>> = witness_evals.par_iter().map(|x| vec![x.c0, x.c1]).collect();
    let mut psi_evals_merkle: Vec<Vec<Fp>> = psi_evals.par_iter().map(|x| vec![x.c0, x.c1]).collect();
    let k: u32 = (((s_ord as f32).log2()).ceil()) as u32;
    witness_evals_merkle =vec![witness_evals_merkle, vec![vec![Fp::from(0),Fp::from(0)]; 2u64.pow(k) as usize - s_ord as usize]].concat();
    psi_evals_merkle =vec![psi_evals_merkle, vec![vec![Fp::from(0),Fp::from(0)]; 2u64.pow(k) as usize - s_ord as usize]].concat();
    
    let w_merkle_slice: Vec<&[Fp]> = witness_evals_merkle.par_iter().map(|x| x.as_slice()).collect();
    let psi_merkle_slice: Vec<&[Fp]> = psi_evals_merkle.par_iter().map(|x| x.as_slice()).collect();
    // Merkle tree of witness evaluations on E
    let witness_mtree: FieldMT = FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        w_merkle_slice,
    )
    .unwrap();
    // Merkle tree of psi evaluations on E
    let psi_mtree: FieldMT = FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        psi_merkle_slice
    )
    .unwrap();
println!("Checkpoint 1");
    let mut roots: Vec<Fp> = vec![witness_mtree.root(), psi_mtree.root()];
    //let roots: Vec<Fp> = vec![Fp::from(0)];
    let alpha_1: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots[..2].to_vec().clone()).unwrap(), Fp::from(0));
    let alpha_2: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_1.c0]).unwrap(), Fp::from(0));
    let alpha_3: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_2.c0]).unwrap(), Fp::from(0));
    let alpha_4: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_3.c0]).unwrap(), Fp::from(0));
    
    // Compute C(x)
    let c1: DensePolynomial<F> = initial_poly(&y_start, b_witness.clone());
    let c2: DensePolynomial<F> = mod_poly_poly(&b_witness, &b_witness_plus, n, g);
    let c3: DensePolynomial<F> = psi_poly(&b_witness, &b_witness_plus_plus, &psi, n, g);
    let c4: DensePolynomial<F> = final_poly(&y_end, b_witness.clone(), g, n as u64);

    let c: DensePolynomial<F> = compute_c(c1, c2, c3, c4, &vec![alpha_1, alpha_2, alpha_3, alpha_4], &n);

    // Evaluate C(x) on E and commit
    let c_evals: Vec<F> = D_0.clone().into_par_iter()
                                     .map(|x| c.evaluate(&x))
                                     .collect();
    let mut c_evals_merkle: Vec<Vec<Fp>> = c_evals.par_iter()
                                            .map(|x| vec![x.c0, x.c1])
                                            .collect();
                                            
    c_evals_merkle =vec![c_evals_merkle, vec![vec![Fp::from(0),Fp::from(0)]; 2u64.pow(k) as usize - s_ord as usize]].concat();
    let c_merkle_slice: Vec<&[Fp]> = c_evals_merkle.par_iter().map(|x| x.as_slice()).collect();
    
    // Merkle tree of witness evaluations on E
    let c_mtree: FieldMT = FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        c_merkle_slice
    )
    .unwrap();
    println!("Checkpoint 2");
    // Compute the evaluations of the constraint polynomial C(x) on E
   
    roots.push(c_mtree.root());

    // Compute the values of the respective polynomials at the challenge z
    let z: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0));

    let gz: F = g*z;
    let ggz: F = g*gz;
    let witness_y: F = b_witness.evaluate(&z);
    let witness_y_plus: F = b_witness.evaluate(&gz);
    let witness_y_plus_plus: F = b_witness.evaluate(&ggz);
    let psi_y: F = psi.evaluate(&z);
    let c_y: F = c.evaluate(&z);
    //let T: u64 = u64::try_from(n).unwrap();
    let E: usize = 32*n/9;

    //assert_eq!(witness_y_plus, b_witness.evaluate(&gz));
    
    let challenge_vals: Vec<F> = vec![witness_y, witness_y_plus, witness_y_plus_plus, psi_y, c_y];
    
    // Finally, create and commit to the composition polynomial P(x)
    let mut zeta_vec: Vec<Fp> = vec![poseidon::CRH::<Fp>::evaluate(&params, vec![z.c0]).unwrap()];
    for _ in 0..4{
        zeta_vec.push(poseidon::CRH::<Fp>::evaluate(&params, zeta_vec.clone()).unwrap());
    }
    let p: DensePolynomial<F> = DensePolynomial { coeffs: vec![vec![F::from(0); E-n], vec![F::new(zeta_vec[0], Fp::from(0))]].concat() }
                                .naive_mul(&(b_witness.clone() + DensePolynomial { coeffs: vec![-challenge_vals[0]] }))
                                .div(&DensePolynomial { coeffs: vec![-z, F::from(1)]})
                                + 
                                DensePolynomial { coeffs: vec![vec![F::from(0); E-n], vec![F::new(zeta_vec[1], Fp::from(0))]].concat() }
                                .naive_mul(&(b_witness.clone() + DensePolynomial { coeffs: vec![-challenge_vals[1]] }))
                                .div(&DensePolynomial { coeffs: vec![-gz, F::from(1)]})
                                + 
                                DensePolynomial { coeffs: vec![vec![F::from(0); E-n], vec![F::new(zeta_vec[2], Fp::from(0))]].concat() }
                                .naive_mul(&(b_witness.clone() + DensePolynomial { coeffs: vec![-challenge_vals[2]] }))
                                .div(&DensePolynomial { coeffs: vec![-ggz, F::from(1)]})
                                +
                                DensePolynomial { coeffs: vec![vec![F::from(0); E-n-1], vec![F::new(zeta_vec[3], Fp::from(0))]].concat() }
                                .naive_mul(&(psi + DensePolynomial { coeffs: vec![-challenge_vals[3]] }))
                                .div(&DensePolynomial { coeffs: vec![-z, F::from(1)]})
                                +
                                DensePolynomial { coeffs:  vec![F::new(zeta_vec[4], Fp::from(0))]}
                                .naive_mul(&(c + DensePolynomial { coeffs: vec![-challenge_vals[4]] }))
                                .div(&DensePolynomial { coeffs: vec![-z, F::from(1)]});


    let (paths_fri, points_fri, roots_fri, indices) = fri_prove(p, l_list, s, r, s_ord, rep_param);
    
println!("Checkpoint 4");
let mut witness_query_vals: Vec<F> = vec![];
let mut witness_plus_query_vals: Vec<F> = vec![];
let mut witness_plus_plus_query_vals: Vec<F> = vec![];
let mut psi_query_vals: Vec<F> = vec![];
let mut c_query_vals: Vec<F> = vec![];

let mut witness_query_path: Vec<FieldPath> = vec![];
let mut witness_plus_query_path: Vec<FieldPath> = vec![];
let mut witness_plus_plus_query_path: Vec<FieldPath> = vec![];
let mut psi_query_path: Vec<FieldPath> = vec![];
let mut c_query_path: Vec<FieldPath> = vec![];
println!("Index 0: {}", indices[0]);
let plus_index: usize = (s_ord as usize)/n;
for index in indices.iter() {
    let z_ind: usize = *index;
    let gz_ind: usize = (z_ind + plus_index) % (s_ord as usize);
    let ggz_ind: usize = (z_ind + 2*plus_index) % (s_ord as usize);
    witness_query_vals.push(witness_evals[z_ind]);
    witness_plus_query_vals.push(witness_evals[gz_ind]);
    witness_plus_plus_query_vals.push(witness_evals[ggz_ind]);
    psi_query_vals.push(psi_evals[z_ind]);
    c_query_vals.push(c_evals[z_ind]);

    witness_query_path.push(witness_mtree.generate_proof(z_ind).unwrap());
    witness_plus_query_path.push(witness_mtree.generate_proof(gz_ind).unwrap());
    witness_plus_plus_query_path.push(witness_mtree.generate_proof(ggz_ind).unwrap());
    psi_query_path.push(psi_mtree.generate_proof(z_ind).unwrap());
    c_query_path.push(c_mtree.generate_proof(z_ind).unwrap());
}
let additional_paths: Vec<Vec<FieldPath>> = vec![witness_query_path, witness_plus_query_path, witness_plus_plus_query_path, psi_query_path, c_query_path];
let additional_points: Vec<Vec<F>> = vec![witness_query_vals, witness_plus_query_vals, witness_plus_plus_query_vals, psi_query_vals, c_query_vals];

(challenge_vals, roots_fri, roots, paths_fri, additional_paths, points_fri, additional_points)


}

pub fn verify(challenges: Vec<F>, roots_fri: Vec<Fp>, roots: Vec<Fp>, paths_fri: Vec<FieldPath>, additional_paths: Vec<Vec<FieldPath>>, points_fri: Vec<F>, additional_points: Vec<Vec<F>>,
            g: F, s: F, r: F, n: &u64, s_ord: u64, y_start: &F, y_end: &F, l_list: Vec<usize>, rep_param: usize) -> bool {
                // Compute z, alphas and zetas 
            let params = poseidon_parameters(); 
            let leaf_crh_params = params.clone();
            let two_to_one_params = params.clone();
            let z: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0));
            let gz: F = g*z;
            let ggz: F = g*gz;
            let alpha_1: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots[..2].to_vec().clone()).unwrap(), Fp::from(0));
            let alpha_2: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_1.c0]).unwrap(), Fp::from(0));
            let alpha_3: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_2.c0]).unwrap(), Fp::from(0));
            let alpha_4: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_3.c0]).unwrap(), Fp::from(0));
            let mut zeta_vec: Vec<Fp> = vec![poseidon::CRH::<Fp>::evaluate(&params, vec![z.c0]).unwrap()];
            
            for _ in 0..4{
                zeta_vec.push(poseidon::CRH::<Fp>::evaluate(&params, zeta_vec.clone()).unwrap());
            }
            // Check that the FRI queries are correct
            let (points_first, indices_first) = fri_verify(paths_fri, points_fri.clone(), roots_fri, l_list.clone(), s.clone(), r.clone(), s_ord.clone(), rep_param as u8);
            // Check that the challenges were computed correctly
            let E: u64 = 32*n/9;
            let c1: F = initial_challenge(y_start, &challenges[0], &z);
            let c2: F = mod_challenge(&challenges[0], &challenges[1], &z, &g, &n);
            let c3: F = psi_challenge(&challenges[0], &challenges[2], &challenges[3], &z, n, &g);
            let c4: F = final_challenge(y_end, &challenges[0], &z, n, &g);
            
            let asserted_c: F = alpha_1*z.pow(&[E-n-2])*c1 
                                + alpha_2*z.pow(&[E-3*n-13])*c2 
                                + alpha_3*z.pow(&[E-n-5])*c3
                                + alpha_4*z.pow(&[E-n-2])*c4;
            assert_eq!(asserted_c, challenges[4]);

            println!("Index 0: {}", indices_first[0]);
            // Check consistency between P(x) in FRI and the committed-to polynomials
            for (i, index) in indices_first.iter().enumerate() {
                let x_0: F = r*s.pow(&[*index as u64]);
                let witness_val: F = additional_points[0][i];
                let witness_plus_val: F = additional_points[1][i];
                let witness_plus_plus_val: F = additional_points[2][i];
                let psi_val: F = additional_points[3][i];
                let c_val: F = additional_points[4][i];
                let asserted_p: F = F::new(zeta_vec[0], Fp::from(0))*x_0.pow(&[E-n])*(witness_val-challenges[0])/(x_0 - z)
                                    + F::new(zeta_vec[1], Fp::from(0))*x_0.pow(&[E-n])*(witness_plus_val-challenges[1])/(x_0 - gz)
                                    + F::new(zeta_vec[2], Fp::from(0))*x_0.pow(&[E-n])*(witness_plus_plus_val-challenges[2])/(x_0 - ggz)
                                    + F::new(zeta_vec[3], Fp::from(0))*x_0.pow(&[E-n-1])*(psi_val-challenges[3])/(x_0 - z)
                                    + F::new(zeta_vec[4], Fp::from(0))*(c_val-challenges[4])/(x_0 - z);
                assert_eq!(asserted_p, points_first[i]);
                // Verify Merkle Paths
                assert!(
                    additional_paths[0][i].verify(
                        &leaf_crh_params,
                        &two_to_one_params,
                        &roots[0],
                        [witness_val.c0, witness_val.c1]
                    ).unwrap());

                assert!(
                    additional_paths[1][i].verify(
                        &leaf_crh_params,
                        &two_to_one_params,
                        &roots[0],
                        [witness_plus_val.c0, witness_plus_val.c1]
                    ).unwrap());

                assert!(
                    additional_paths[2][i].verify(
                        &leaf_crh_params,
                        &two_to_one_params,
                        &roots[0],
                        [witness_plus_plus_val.c0, witness_plus_plus_val.c1]
                    ).unwrap());

                assert!(
                    additional_paths[3][i].verify(
                        &leaf_crh_params,
                        &two_to_one_params,
                        &roots[1],
                        [psi_val.c0, psi_val.c1]
                    ).unwrap());
                assert!(
                    additional_paths[4][i].verify(
                        &leaf_crh_params,
                        &two_to_one_params,
                        &roots[2],
                        [c_val.c0, c_val.c1]
                    ).unwrap());
            }
            true}
pub fn mod_challenge(x: &F, y: &F, z: &F, g: &F, T: &u64) -> F {
    let eval: F = x*x*x+y*y*y-x*x*y*y+F::from(1488u128)
    *(x*x*y+y*y*x)-F::from(162000u128)*
    (x*x+y*y)+F::from(40773375u128)*x*y
    +F::from(8748000000u128)*(x+y)-
    F::from(157464000000000u128);

    (z-g.pow(&[*T-1]))*eval / (z.pow(&[*T])-F::from(1))
}

//Returns (p(x)-y_0)/(x - 1)
    pub fn initial_poly(y_0: &F, p: DensePolynomial<F>) -> DensePolynomial<F> {
        (p + DensePolynomial{coeffs: vec![-*y_0]}).div(&DensePolynomial{coeffs: vec![F::from(-1), F::from(1)]})
    }
    //Returns (p(x)-y_end)/(x - g^(T-1))
    pub fn final_poly(y_end: &F, p: DensePolynomial<F>, g: F, T: u64) -> DensePolynomial<F> {
        (p + DensePolynomial{coeffs: vec![-*y_end]}).div(&DensePolynomial{coeffs: vec![-g.pow(&[T-1]), F::from(1)]})
    }

    //Returns the composite polynomial
    pub fn compute_c(c1: DensePolynomial<F>, c2: DensePolynomial<F>, c3: DensePolynomial<F>, c4: DensePolynomial<F>, alphas: &Vec<F>, T: &usize) -> DensePolynomial<F> {
        let E: usize = 32*T/9;
        let deg_1: usize = E-T-2;
        let deg_2: usize = E-T-5;
        let deg_3: usize = E-3*T - 13;

        c1.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_1], vec![alphas[0]]].concat()})
        + c2.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_3], vec![alphas[1]]].concat()})
        + c3.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_2], vec![alphas[2]]].concat()})
        + c4.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_1], vec![alphas[3]]].concat()})
    }
    // Returns (x-g^(T-1))*Phi_2(p(x), q(x))/(x^n - 1)
   pub fn mod_poly_poly(p: &DensePolynomial<F>, q: &DensePolynomial<F>, T: usize, g: F) -> DensePolynomial<F> {
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

        temp.naive_mul(&DensePolynomial{coeffs: vec![-g.pow(&[T as u64-1]), F::from(1)]}).div(&DensePolynomial{ coeffs: [vec![-F::from(1)], vec![F::from(0); T-1], vec![F::from(1)]].concat()})
    }
    
    pub fn initial_challenge(y_0: &F, eval: &F, x_0: &F) -> F {

        (eval - y_0) / (x_0-F::from(1))

    }
    pub fn final_challenge(y_end: &F, eval: &F, x_0: &F, n: &u64, g: &F) -> F {

        (eval - y_end) / (x_0-g.pow(&[*n-1]))
    }
    //Returns ((x-g^(T-2))*(x-g^(T-1))*(p(x)-q(x))*psi(x)-1)/(x^T-1)
   pub fn psi_poly(p: &DensePolynomial<F>, q: &DensePolynomial<F>, psi: &DensePolynomial<F>, T: usize, g: F) -> DensePolynomial<F> {
        let diff:DensePolynomial<F> = p.sub(q);
        let g_pow: F =  g.pow(&[T as u64-2]);
        let x_1_poly: DensePolynomial<F> = DensePolynomial{coeffs: vec![-g_pow, F::from(1)]}.naive_mul(&DensePolynomial{coeffs: vec![-g_pow*g, F::from(1)]});
        
        x_1_poly.naive_mul(&(diff.naive_mul(psi)+DensePolynomial{coeffs:vec![F::from(-1)]})).div(&DensePolynomial{ coeffs: [vec![-F::from(1)], vec![F::from(0); T-1], vec![F::from(1)]].concat()})

        }
    
    pub fn psi_challenge(y_witness: &F, y_witness_plusplus: &F, y_psi: &F, x_0: &F, n: &u64, g: &F) -> F {
        let g_pow: F =  g.pow(&[*n as u64-2]);
        let g_prefactor: F = (x_0 - g_pow) * (x_0 - g_pow*g);
        
        g_prefactor*((y_witness-y_witness_plusplus)*y_psi-F::from(1))/(x_0.pow(&[*n])-F::from(1))
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