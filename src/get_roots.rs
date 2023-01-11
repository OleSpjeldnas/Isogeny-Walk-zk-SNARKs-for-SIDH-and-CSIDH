use super::*;

// Computes a valid 2-isogeny walk of length n
pub fn mod_poly_eval_at_y_div_by_x_minus_z(y: F, z: F) -> Vec<F> {
    let mod_poly_at_y: Vec<F> = vec![
        -F::from(157464000000000u128) + F::from(8748000000u128) * y - F::from(162000u128) * y.square() + y.pow([3]),
        F::from(8748000000u128) + F::from(40773375u128) * y + F::from(1488u128) * y.square(),
        F::from(1488u128) * y - F::from(162000u128) - y.square(),
        F::from(1),
    ];
    vec![mod_poly_at_y[1] + z * (mod_poly_at_y[2] + z), mod_poly_at_y[2] + z, F::from(1)]
}

pub fn find_roots(roots_arr: &mut Vec<F>, n: usize) -> Vec<F> {
    for _ in 0..n {
        let l = roots_arr.len();
        let z: F = roots_arr[l - 2];
        let y: F = roots_arr[l - 1];

        let poly: Vec<F> = mod_poly_eval_at_y_div_by_x_minus_z(y, z);
        //let deriv: Vec<F> = vec![poly[1], poly[2]*F::from(2)];
        //roots_arr.push(newton_method(poly, deriv, y));
        //let root: F = -zassenhaus(poly)[0][0];
        //roots_arr.push(root);
        roots_arr.push((-poly[1] + (poly[1].square() - F::from(4u128) * poly[0]).sqrt().unwrap()) / F::from(2u128));
    }
    roots_arr.to_vec()
}
