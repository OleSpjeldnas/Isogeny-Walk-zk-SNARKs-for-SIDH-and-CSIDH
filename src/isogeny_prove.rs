use std::ops::{Div, Mul, Sub};
use rayon::prelude::*;

use super::*;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::{UniformRand, FftField};
use ark_std::test_rng;
use merkle::{FieldMT, poseidon_parameters};
// Witness is the witness polynomial, psi the inverse of w(x)-w(g^2x), g the generator of the interpolation domain, 
//the evaluation domain is E = r<s>. Finally, s_ord is the size of E.
pub fn prove(witness: DensePolynomial<F>, witness_plus: DensePolynomial<F>, witness_plus_plus: DensePolynomial<F>,                 // Challenges, roots, roots_fri, paths_fri, paths, points_fri, points
    psi: DensePolynomial<F>, g: F, s: F, r: F, s_ord: u64, y_start: &F, y_end: &F, l_list: Vec<usize>, rep_param: usize) -> (Vec<F>, Vec<Fp>, Vec<Fp>, Vec<FieldPath>, Vec<Vec<FieldPath>>, Vec<F>, Vec<Vec<F>>){
    let n: usize = witness.coeffs.len();
    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    //Blind the witness
    let b_witness: DensePolynomial<F> =  witness + DensePolynomial { coeffs: vec![vec![a, b, c], vec![F::from(0); n-1], vec![a, b, c]].concat()};
    let b_witness_plus: DensePolynomial<F> =  witness_plus + DensePolynomial { coeffs: vec![vec![a, b, c], vec![F::from(0); n-1], vec![a, b, c]].concat()};
    let b_witness_plus_plus: DensePolynomial<F> =  witness_plus_plus + DensePolynomial { coeffs: vec![vec![a, b, c], vec![F::from(0); n-1], vec![a, b, c]].concat()};

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
    let alpha_1: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0));
    let alpha_2: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_1.c0]).unwrap(), Fp::from(0));
    let alpha_3: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_2.c0]).unwrap(), Fp::from(0));
    let alpha_4: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_3.c0]).unwrap(), Fp::from(0));

    // Compute C(x)
    let c1: DensePolynomial<F> = initial_poly(&y_start, b_witness.clone());
    let c4: DensePolynomial<F> = final_poly(&y_end, b_witness.clone(), g, n as u64);
    let c2: DensePolynomial<F> = mod_poly_poly(&b_witness, &b_witness_plus, n, g);
    let c3: DensePolynomial<F> = psi_poly(&b_witness, &b_witness_plus_plus, &psi, n, g);

    let c: DensePolynomial<F> = compute_c(c1, c2, c3, c4, &vec![alpha_1, alpha_2, alpha_3, alpha_4], &n);

    // Evaluate C(x) on E and commit
    let c_evals: Vec<F> = D_0.clone().into_par_iter()
                                     .map(|x| c.evaluate(&x))
                                     .collect();
    let mut c_evals_merkle: Vec<Vec<Fp>> = c_evals.par_iter()
                                            .map(|x| vec![x.c0, x.c1])
                                            .collect();
    let k: u32 = (((s_ord as f32).log2()).ceil()) as u32;
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
    roots.push(c_mtree.root());
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
    let c_y: F = F::from(0);
    //let T: u64 = u64::try_from(n).unwrap();
    let E: usize = 32*n/9;
    let deg_1: u64 = (E-n-2).try_into().unwrap();
    let deg_2: u64 = (E-n-5).try_into().unwrap();

    //let init_term: F = (witness_y - y_start) / (z.pow(&[T])-F::from(1));
    //let z_pow_T_minus_1: F = z.pow(&[T]) - F::from(1);
    //let g_T_2: F = g.pow(&[T-2]);
    //let g_T_1: F = g*g_T_2;
    // C2(z)
    //let mod_term: F = (z-g_T_1)*mod_poly(witness_y, witness_y_plus) / z_pow_T_minus_1;
    // C3(z)
    //let psi_term: F = (z-g_T_2)*(z-g_T_1)* ((witness_y - witness_y_plus_plus)*psi_y - F::from(1)) / z_pow_T_minus_1;
    //let c_y: F = alpha_1*z.pow(&[deg_1])*init_term + mod_term + alpha_2*z.pow(&[deg_2])*psi_term;

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
                                .div(&DensePolynomial { coeffs: vec![-z, F::from(1)]});
                                //+
                                //DensePolynomial { coeffs:  vec![F::new(zeta_vec[4], Fp::from(0))]}
                                //.naive_mul(&(c + DensePolynomial { coeffs: vec![-challenge_vals[4]] }))
                                //.div(&DensePolynomial { coeffs: vec![-z, F::from(1)]});


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

let plus_index: usize = (s_ord as usize)/n;
for index in indices.iter() {
    witness_query_vals.push(witness_evals[*index]);
    witness_plus_query_vals.push(witness_evals[(*index + plus_index) % n]);
    witness_plus_plus_query_vals.push(witness_evals[(*index + 2*plus_index) % n]);
    psi_query_vals.push(psi_evals[*index]);
    c_query_vals.push(c_evals[*index]);

    witness_query_path.push(witness_mtree.generate_proof(*index).unwrap());
    witness_plus_query_path.push(witness_mtree.generate_proof((*index + plus_index) % n).unwrap());
    witness_plus_plus_query_path.push(witness_mtree.generate_proof((*index + 2*plus_index) % n).unwrap());
    psi_query_path.push(psi_mtree.generate_proof(*index).unwrap());
    c_query_path.push(c_mtree.generate_proof(*index).unwrap());
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
            let c1: F = initial_challenge(y_start, &challenges[0], &z, n);
            let c2: F = final_challenge(y_end, &challenges[0], &z, n, &g);
            let c3: F = mod_challenge(&challenges[0], &challenges[1], &z, &g, &n);
            let c4: F = psi_challenge(&challenges[0], &challenges[2], &challenges[3], &z, n, &g);
            
            let asserted_c: F = alpha_1*z.pow(&[E-n-2])*c1 
                                + alpha_3*z.pow(&[E-n-2])*c2 
                                + alpha_4*z.pow(&[E-3*n-13])*c3
                                + alpha_2*z.pow(&[E-n-5])*c4;
            assert_eq!(asserted_c, challenges[4]);

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
    fn mod_challenge(x: &F, y: &F, z: &F, g: &F, T: &u64) -> F {
        let eval: F = x*x*x+y*y*y-x*x*y*y+F::from(1488u128)
        *(x*x*y+y*y*x)-F::from(162000u128)*
        (x*x+y*y)+F::from(40773375u128)*x*y
        +F::from(8748000000u128)*(x+y)-
        F::from(157464000000000u128);

        (z-g.pow(&[*T-1]))*eval / (z.pow(&[*T])-F::from(1))
    }
    //evaluates Phi_2(witness(x), witness_plus(x))
    fn mod_poly_const_eval_at_x(x: &F, witness: &DensePolynomial<F>, g: &F) -> F {
        let z: F = witness.evaluate(&x);
        let y: F = witness.evaluate(&(g*x));
        let z_squared: F = z*z;
        let y_squared: F = y*y;
        let z_cubed: F = z_squared*z;
        let y_cubed: F = y_squared*y;
        let z_squared_y: F = z_squared*y;
        let y_squared_z: F = y_squared*z;
        let yz: F = y*z;

        z_cubed+y_cubed-z_squared*y_squared+F::from(1488u128)
        *(z_squared_y+y_squared_z)-F::from(162000u128)*
        (z_squared+y_squared)+F::from(40773375u128)*yz
        +F::from(8748000000u128)*(z+y)-
        F::from(157464000000000u128)
    }

    //evaluates (witness(x) - y_0) / (x^n - 1) at some point x
    fn initial_constraint(x: &F, y_0: &F, witness: &DensePolynomial<F>, n: &u64) -> F {

        (witness.evaluate(&x) - y_0) / (x.pow(&[*n])-F::from(1))
    }
    fn psi_constraint(x: &F, witness: &DensePolynomial<F>, psi: &DensePolynomial<F>, g: &F, n: &u64) -> F {
        let z: F = witness.evaluate(&x);
        let y: F = witness.evaluate(&(g.square()*x));
        let w: F = psi.evaluate(&x);
        let g_pow: F =  g.pow(&[*n-2]);
        let x_1: F = x-g_pow;
        let x_2: F = x_1-g_pow*g;

        (x_1*x_2*(z-y)*w-F::from(1)) / (x.pow(&[*n])-F::from(1))

    }
    
    fn compute_c_evaluations(witness_evals: &Vec<F>, psi_evals: &Vec<F>, n: &usize, s_ord: &usize, r: &F, s: &F, y_0: &F, g: &F, alphas: &Vec<F>) -> Vec<Vec<Fp>>{
        let mut c_evals: Vec<Vec<Fp>> = Vec::new();
        let T: u64 = u64::try_from(*n).unwrap();
        let E: u64 = 3*T + 13;
        let deg_1: u64 = E-T-2;
        let deg_2: u64 = E-T-5;
        for i in 0..*s_ord {
            // Set x := r*s^i
            let x: F = r*s.pow(&[i as u64]);
            // C1(x)
            let init_term: F = (witness_evals[i] - y_0) / (x.pow(&[*n as u64])-F::from(1));
            let x_pow_T_minus_1: F = x.pow(&[*n as u64]) - F::from(1);
            let g_T_2: F = g.pow(&[*n as u64-2]);
            let g_T_1: F = g*g_T_2;
            // C2(x)
            let mod_term: F = (x-g_T_1)*mod_poly(witness_evals[i], witness_evals[i+1]) / x_pow_T_minus_1;
            // C3(x)
            let psi_term: F = (x-g_T_2)*(x-g_T_1)* ((witness_evals[i] - witness_evals[i+2])*psi_evals[i] - F::from(1)) / x_pow_T_minus_1;
            let C_term: F = alphas[0]*x.pow(&[deg_1])*init_term + mod_term + alphas[1]*x.pow(&[deg_2])*psi_term;

            c_evals.push(vec![C_term.c0, C_term.c1]);
        }
        let k = (((*s_ord as f32).log2()).ceil()) as u8;
    for _ in 0..(2u64.pow(k.into())-u64::try_from(*s_ord).unwrap()) {
        c_evals.push(vec![Fp::from(0),Fp::from(0)]);
    }
        c_evals
    }

    //Returns (p(x)-y_0)/(x - 1)
    fn initial_poly(y_0: &F, p: DensePolynomial<F>) -> DensePolynomial<F> {
        (p.clone() + DensePolynomial{coeffs: vec![-*y_0]}).div(&DensePolynomial{coeffs: vec![F::from(-1), F::from(1)]})
    }
    //Returns (p(x)-y_end)/(x - g^(T-1))
    fn final_poly(y_end: &F, p: DensePolynomial<F>, g: F, T: u64) -> DensePolynomial<F> {
        (p + DensePolynomial{coeffs: vec![-*y_end]}).div(&DensePolynomial{coeffs: vec![-g.pow(&[T-1]), F::from(1)]})
    }

    //Returns the composite polynomial
    fn compute_c(c1: DensePolynomial<F>, c2: DensePolynomial<F>, c3: DensePolynomial<F>, c4: DensePolynomial<F>, alphas: &Vec<F>, T: &usize) -> DensePolynomial<F> {
        let E: usize = 32*T/9;
        let deg_1: usize = E-T-2;
        let deg_2: usize = E-T-5;
        let deg_3: usize = 3*T + 13;

        c1.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_1], vec![alphas[0]]].concat()})
        + c2.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_3], vec![alphas[3]]].concat()})
        + c3.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_2], vec![alphas[1]]].concat()})
        + c4.naive_mul(&DensePolynomial{coeffs: vec![vec![F::from(0); deg_1], vec![alphas[2]]].concat()})
    }
    // Returns (x-g^(T-1))*Phi_2(p(x), q(x))/(x^n - 1)
    fn mod_poly_poly(p: &DensePolynomial<F>, q: &DensePolynomial<F>, T: usize, g: F) -> DensePolynomial<F> {
        let p_squared: DensePolynomial<F> = p.naive_mul(p);
        let q_squared: DensePolynomial<F> = q.naive_mul(q);
        let p_cubed: DensePolynomial<F> = p_squared.naive_mul(p);
        let q_cubed: DensePolynomial<F> = q_squared.naive_mul(q);
        let p_squared_q: DensePolynomial<F> = p_squared.naive_mul(q);
        let q_squared_p: DensePolynomial<F> = q_squared.naive_mul(p);
        let pq: DensePolynomial<F> = p.naive_mul(q);

        let temp: DensePolynomial<F> = p_cubed+q_cubed.sub(&p_squared_q.naive_mul(q))+DensePolynomial{coeffs: vec![F::from(1488u128)]}.naive_mul(&(p_squared_q+q_squared_p))
        +DensePolynomial{coeffs: vec![F::from(-162000i32)]}.naive_mul(&(p_squared+q_squared))
        +DensePolynomial{coeffs: vec![F::from(40773375u128)]}.naive_mul(&pq)
        +DensePolynomial{coeffs: vec![F::from(8748000000u128)]}.naive_mul(&(p+q))+
        DensePolynomial{coeffs: vec![F::from(-157464000000000i64)]};

        temp.naive_mul(&DensePolynomial{coeffs: vec![-g.pow(&[T as u64-1]), F::from(1)]}).div(&DensePolynomial{ coeffs: [vec![F::from(1)], vec![F::from(0); T-1], vec![F::from(1)]].concat()})
    }
    
    fn initial_challenge(y_0: &F, eval: &F, x_0: &F, n: &u64) -> F {

        (eval - y_0) / (x_0.pow(&[*n])-F::from(1))

    }
    fn final_challenge(y_eng: &F, eval: &F, x_0: &F, n: &u64, g: &F) -> F {

        (eval - y_eng) / (x_0.pow(&[*n])-g.pow(&[*n-1]))
    }
    //Returns ((x-g^(T-2))*(x-g^(T-1))*(p(x)-q(x))*psi(x)-1)/(x^T-1)
    fn psi_poly(p: &DensePolynomial<F>, q: &DensePolynomial<F>, psi: &DensePolynomial<F>, T: usize, g: F) -> DensePolynomial<F> {
        let diff:DensePolynomial<F> = p.sub(q);
        let g_pow: F =  g.pow(&[T as u64-2]);
        let x_1_poly: DensePolynomial<F> = DensePolynomial{coeffs: vec![-g_pow, F::from(1)]}.naive_mul(&DensePolynomial{coeffs: vec![-g_pow*g, F::from(1)]});
        
        x_1_poly.naive_mul(&(diff.naive_mul(psi)+DensePolynomial{coeffs:vec![F::from(-1)]})).div(&DensePolynomial{ coeffs: [vec![F::from(1)], vec![F::from(0); T-1], vec![F::from(1)]].concat()})

        }
    
    fn psi_challenge(y_witness: &F, y_witness_plusplus: &F, y_psi: &F, x_0: &F, n: &u64, g: &F) -> F {
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