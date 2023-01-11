use super::*;

pub fn mod_poly_eval_at_y_div_by_x_minus_z(y: F, z: F) -> Vec<F> {
    let mod_poly_at_y: Vec<F> = vec![
        -F::from(157464000000000u128) + F::from(8748000000u128) * y - F::from(162000u128) * y.square() + y.pow(&[3]),
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

fn zassenhaus(poly: Vec<F>) -> Vec<Vec<F>> {
    let n = poly.len() - 1;
    let mut factors = Vec::new();

    if n == 1 {
        // The polynomial is linear, so it has only one factor (x - poly[0]).
        factors.push(vec![-poly[0], F::default()]);
        return factors;
    }

    let mut g = vec![F::default(); n + 1];
    g[n] = -poly[0] / poly[n];
    let mut h = vec![F::default(); n];
    h[n - 1] = F::default();

    for i in (1..n).rev() {
        let mut s = poly[i];
        for j in (i + 1)..(n + 1) {
            s += g[j] * poly[j - i - 1];
        }
        g[i] = -s / poly[n];
        for j in (0..i).rev() {
            h[j] = h[j] - g[i] * h[i - j - 1];
        }
    }

    factors.push(h);
    factors.push(g);

    factors
}

fn newton_method(f: Vec<F>, df: Vec<F>, x0: F) -> F {
    let mut x: F = x0;
    loop {
        let fx: F = f[0] + f[1] * x + f[2] * x.square();
        let dfx: F = df[0] + df[1] * x;
        let x_new: F = x - fx / dfx;
        if (x - x_new) < F::from(1) || -(x - x_new) < F::from(1) {
            return x_new;
        }
        x = x_new;
    }
}
