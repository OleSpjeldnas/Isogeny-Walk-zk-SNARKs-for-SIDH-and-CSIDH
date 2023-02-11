pub fn isogeny_run () -> () {
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
}