use super::*;
use ark_poly::polynomial::univariate::DensePolynomial;
use merkle::{FieldMT, poseidon_parameters};
//use core::slice::SlicePattern;
use std::ops::{Div};
// Witness is the witness polynomial, psi the inverse of w(x)-w(g^2x), g the generator of the interpolation domain, t the length of each walk
//the evaluation domain is E = r<s>. Finally, s_ord is the size of E.   ->    (Challenges, roots, roots_fri, paths_fri, paths, points_fri, points)
pub fn prove(mod_polys: &Vec<Vec<(u32, Vec<(u32, K)>)>>, witnesses: Vec<DensePolynomial<K>>, witnesses_g: Vec<DensePolynomial<K>>, g: K, t: &u8, s: K, r: K, s_ord: &u64, y_start: &K, y_end: &K, l_list: &Vec<usize>, rep_param: usize, grinding_param: u8)
 -> (Vec<K>, Vec<Kp>, Vec<Kp>, Vec<FieldPath>, Vec<Vec<FieldPath>>, Vec<K>, Vec<Vec<K>>)
{
    let main_witness: DensePolynomial<K> = witnesses[1].clone();
    let T: usize = (t*74) as usize;

    let params = poseidon_parameters();
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();


    let D_0: Vec<K> = (0..*s_ord).into_iter()
                                .map(|i| r*s.pow([i]))
                                .collect();
    let witness_evals: Vec<K> = D_0.clone().into_iter()
                                   .map(|x| main_witness.evaluate(&x))
                                   .collect();
    let mut witness_evals_merkle: Vec<Vec<Kp>> = witness_evals.
    iter()
    .map(|x| vec![x.c0, x.c1])
    .collect();
    
    let k: u32 = (((*s_ord as f32).log2()).ceil()) as u32;
    witness_evals_merkle =vec![witness_evals_merkle, vec![vec![Kp::from(0),Kp::from(0)]; 2u64.pow(k) as usize - *s_ord as usize]].concat();
    let merkle_slice: Vec<&[Kp]> = witness_evals_merkle.
                                    iter()
                                    .map(|x| x.as_slice())
                                    .collect();
    let witness_mtree: FieldMT = FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        merkle_slice,
    )
    .unwrap();
    let mut roots: Vec<Kp> = vec![witness_mtree.root()];
    let alpha_1: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, roots.clone()).unwrap(), Kp::from(0));
    let alpha_2: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, vec![alpha_1.c0]).unwrap(), Kp::from(0));
    let alpha_3: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, vec![alpha_2.c0]).unwrap(), Kp::from(0));
    // Compute constraint polynomials
    let c1: DensePolynomial<K> = initial_poly_constraint(&y_start, &main_witness);
    let c2: DensePolynomial<K> = final_poly_constraint(&y_end, &main_witness, &g, &T);
    let c3: DensePolynomial<K> = transition_constraint(mod_polys, &witnesses, &witnesses_g, &g, &(*t as usize));

    // 2 x 3 x 5 x 7 x 11 x 13 x 17
    //let E: usize = 510510;
    let E: usize = T*12;
    let deg_1: usize = T-1;
    let deg_2: usize = T-1;
    //let deg_3: usize = T*(1+2*587)-(*t as usize);
    let deg_3: usize = T*(1+2*5)-(*t as usize);
    // Compute composite constraint polynomial C(x)
    //let c: DensePolynomial<K> = c1.naive_mul(&DensePolynomial{coeffs: vec![vec![K::from(0); E-deg_1], vec![alpha_1]].concat()})
    //+ c2.naive_mul(&DensePolynomial{coeffs: vec![vec![K::from(0); E-deg_2], vec![alpha_2]].concat()})
    //+ c3.naive_mul(&DensePolynomial{coeffs: vec![vec![K::from(0); E-deg_3], vec![alpha_3]].concat()});

    let c: DensePolynomial<K> = c3.naive_mul(&DensePolynomial{coeffs: vec![vec![K::from(0); E-deg_3], vec![alpha_3]].concat()});
    println!("Checkpoint 0");

    let c_evals: Vec<K> = D_0.clone().into_iter()
                                   .map(|x| c.evaluate(&x))
                                   .collect();
    let mut c_evals_merkle: Vec<Vec<Kp>> = c_evals.
    iter()
    .map(|x| vec![x.c0, x.c1])
    .collect();
    println!("Checkpoint 1");
    c_evals_merkle =vec![c_evals_merkle, vec![vec![Kp::from(0),Kp::from(0)]; 2u64.pow(k) as usize - *s_ord as usize]].concat();
    let c_merkle_slice: Vec<&[Kp]> = c_evals_merkle.
                                    iter()
                                    .map(|x| x.as_slice())
                                    .collect();
    let c_mtree: FieldMT = FieldMT::new(
        &leaf_crh_params,
        &two_to_one_params,
        c_merkle_slice,
    )
    .unwrap();
    roots.push(c_mtree.root());
    println!("Checkpoint 2");

    // Now use Fiat-SHamir to compute the challenge z
    let z: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, roots.clone()).unwrap(), Kp::from(0));
    let gz: K = g*z;
    let witness_y: K = main_witness.evaluate(&z);
    let witness_y_plus: K = main_witness.evaluate(&gz);
    let c_y: K = c.evaluate(&z);
    let challenge_vals: Vec<K> = vec![witness_y, witness_y_plus, c_y];
    // Compute P(x)
    let mut zeta_vec: Vec<Kp> = vec![poseidon::CRH::<Kp>::evaluate(&params, vec![z.c0]).unwrap()];
    for _ in 0..2{
        zeta_vec.push(poseidon::CRH::<Kp>::evaluate(&params, zeta_vec.clone()).unwrap());
    }
    println!("Checkpoint 3");
    let p: DensePolynomial<K> = DensePolynomial { coeffs: vec![vec![K::from(0); E-T], vec![K::new(zeta_vec[0], Kp::from(0))]].concat() }
                                .naive_mul(&(main_witness.clone() + DensePolynomial { coeffs: vec![-challenge_vals[0]] }))
                                .div(&DensePolynomial { coeffs: vec![-z, K::from(1)]})
                                + 
                                DensePolynomial { coeffs: vec![vec![K::from(0); E-T], vec![K::new(zeta_vec[1], Kp::from(0))]].concat() }
                                .naive_mul(&(main_witness + DensePolynomial { coeffs: vec![-challenge_vals[1]] }))
                                .div(&DensePolynomial { coeffs: vec![-gz, K::from(1)]})
                                + 
                                DensePolynomial { coeffs:  vec![K::new(zeta_vec[2], Kp::from(0))]}
                                .naive_mul(&(c + DensePolynomial { coeffs: vec![-challenge_vals[2]] }))
                                .div(&DensePolynomial { coeffs: vec![-z, K::from(1)]});

    let (paths_fri, points_fri, roots_fri, indices) = fri_prove(p.clone(), l_list.to_vec(), s, r, *s_ord, rep_param, grinding_param);

    let mut witness_query_vals: Vec<K> = vec![];
    let mut c_query_vals: Vec<K> = vec![];
    
    let mut witness_query_path: Vec<FieldPath> = vec![];
    let mut c_query_path: Vec<FieldPath> = vec![];
    for index in indices.iter() {
        let z_ind: usize = *index;
        witness_query_vals.push(witness_evals[z_ind]);
        c_query_vals.push(c_evals[z_ind]);
    
        witness_query_path.push(witness_mtree.generate_proof(z_ind).unwrap());
        c_query_path.push(c_mtree.generate_proof(z_ind).unwrap());
                                }
let additional_paths: Vec<Vec<FieldPath>> = vec![witness_query_path, c_query_path];
let additional_points: Vec<Vec<K>> = vec![witness_query_vals, c_query_vals];
println!("Checkpoint 4");

(challenge_vals, roots_fri, roots, paths_fri, additional_paths, points_fri, additional_points)
}

pub fn verify(mod_polys: &Vec<Vec<(u32, Vec<(u32, K)>)>>, challenges: Vec<K>, roots_fri: Vec<Kp>, roots: Vec<Kp>, paths_fri: Vec<FieldPath>, additional_paths: Vec<Vec<FieldPath>>, points_fri: Vec<K>, additional_points: Vec<Vec<K>>,
    g: K, s: K, r: K, s_ord: u64, t: u8, y_start: &K, y_end: &K, l_list: Vec<usize>, rep_param: usize, grinding_param: u8) 
    -> bool 
    {
        // Compute z, alphas and zetas 
    let T: u64  = 74*t as u64;
    let params = poseidon_parameters(); 
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();
    let z: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, roots.clone()).unwrap(), Kp::from(0));
    let gz: K = g*z;
    let alpha_1: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, vec![roots[0]]).unwrap(), Kp::from(0));
    let alpha_2: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, vec![alpha_1.c0]).unwrap(), Kp::from(0));
    let alpha_3: K = K::new(poseidon::CRH::<Kp>::evaluate(&params, vec![alpha_2.c0]).unwrap(), Kp::from(0));

    let mut zeta_vec: Vec<Kp> = vec![poseidon::CRH::<Kp>::evaluate(&params, vec![z.c0]).unwrap()];
            
    for _ in 0..2{
        zeta_vec.push(poseidon::CRH::<Kp>::evaluate(&params, zeta_vec.clone()).unwrap());
    }

    let (points_first, indices_first) = fri_verify(paths_fri, points_fri.clone(), roots_fri, l_list.clone(), s.clone(), r.clone(), s_ord.clone(), rep_param as u8, grinding_param);
    // Check that the challenges were computed correctly
    //let E: u64 = 510510;
    let E: u64 = 74*5*12;
    let deg_1: u64 = T -1;
    let deg_2: u64 = T-1;
    //let deg_3: u64 = T*(1+2*587)-(t as u64);
    let deg_3: u64 = T*(1+2*5)-(t as u64);

    let c1: K = initial_verifier(y_start, &challenges[0], &z);
    let c2: K = final_verifier(&y_end, &challenges[0], &z, &g, &T);
    let c3: K = transition_verifier(mod_polys, &challenges[0], &challenges[1], &z, &g, &(t as usize));

    //let asserted_c: K = alpha_1*z.pow(&[E-deg_1])*c1 
    //                                        + alpha_2*z.pow(&[E-deg_2])*c2 
    //                                        + alpha_3*z.pow(&[E-deg_3])*c3;
    let asserted_c: K =alpha_3*z.pow(&[E-deg_3])*c3;
    
    assert_eq!(asserted_c, challenges[2]);

   for (i, index) in indices_first.iter().enumerate() {
                let x_0: K = r*s.pow(&[*index as u64]);
                
                let witness_val: K = additional_points[0][i];
                let c_val: K = additional_points[1][i];
                let asserted_p: K = K::new(zeta_vec[0], Kp::from(0))*x_0.pow(&[E-T])*(witness_val-challenges[0])/(x_0 - z)
                                    + K::new(zeta_vec[1], Kp::from(0))*x_0.pow(&[E-T])*(witness_val-challenges[1])/(x_0 - gz)
                                    + K::new(zeta_vec[2], Kp::from(0))*(c_val-challenges[2])/(x_0 - z);
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
                        &roots[1],
                        [c_val.c0, c_val.c1]
                    ).unwrap());
            }
            
            true}


fn initial_poly_constraint(y_start: &K, p: &DensePolynomial<K>) -> DensePolynomial<K> {
    (p + &DensePolynomial{coeffs: vec![-*y_start]}).div(&DensePolynomial{coeffs: vec![-K::from(1), K::from(1)]})
}
fn initial_verifier(y_start: &K, eval: &K, x_0: &K) -> K {
    (eval - y_start) / (x_0-K::from(1))
}
//Returns (p(x)-y_end)/(x - g^(T-1))
fn final_poly_constraint(y_end: &K, p: &DensePolynomial<K>, g: &K, T: &usize) -> DensePolynomial<K> {
    (p + &DensePolynomial{coeffs: vec![-*y_end]}).div(&DensePolynomial{coeffs: vec![-g.pow(&[*T as u64 -1]), K::from(1)]})
}
fn final_verifier(y_end: &K, eval: &K, x_0: &K, g: &K, T: &u64) -> K {
    (eval - y_end) / (x_0-g.pow(&[*T-1]))
}
fn transition_constraint(mod_polys: &Vec<Vec<(u32, Vec<(u32, K)>)>>, witnesses: &Vec<DensePolynomial<K>>, witnesses_g: &Vec<DensePolynomial<K>>, g: &K, t: &usize)
 -> DensePolynomial<K> 
{   let mut g_pow = K::from(1);
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

fn transition_verifier(mod_polys: &Vec<Vec<(u32, Vec<(u32, K)>)>>, y: &K, y_g: &K, z: &K, g: &K, t: &usize)
 -> K 
{   let mut ys: Vec<K> = vec![K::from(1)];
    let mut y_gs: Vec<K> = vec![K::from(1)];
    for _ in 0..588 {
        ys.push(ys.last().unwrap() * y);
        y_gs.push(y_gs.last().unwrap() * y_g);
    }
    let mut g_pow = K::from(1);
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