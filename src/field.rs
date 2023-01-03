use ark_ff::{Fp448, MontBackend, Fp2, MontFp, Fp2Config, Field};


#[derive(ark_ff::fp::MontConfig)]
#[modulus = "46125640644791984503063093131192152308860971046482148971221392549658928357566369654107368811182716603910699644514503519777239072769"]
#[generator = "17"]

pub struct FqConfig;
pub type F = Fp448<MontBackend<FqConfig, 7>>;
pub type Fq2 = Fp2<Fq2Config>;

pub struct Fq2Config;

impl Fp2Config for Fq2Config {
    type Fp = F;

    /// NONRESIDUE = -1
    const NONRESIDUE: F = MontFp!("17");

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP2_C1: &'static [F] = &[
        MontFp!("46125640644791984503063093131192152308860971046482148971221392549658928357566369654107368811182716603910699644514503519777239072768"),
        MontFp!("46125640644791984503063093131192152308860971046482148971221392549658928357566369654107368811182716603910699644514503519777239072768")
    ];

    #[inline(always)]
    fn mul_fp_by_nonresidue_in_place(fp: &mut Self::Fp) -> &mut Self::Fp {
        fp.neg_in_place()
    }

    #[inline(always)]
    fn sub_and_mul_fp_by_nonresidue(y: &mut Self::Fp, x: &Self::Fp) {
        *y += x;
    }

    #[inline(always)]
    fn mul_fp_by_nonresidue_plus_one_and_add(y: &mut Self::Fp, x: &Self::Fp) {
        *y = *x;
    }

    fn mul_fp_by_nonresidue_and_add(y: &mut Self::Fp, x: &Self::Fp) {
        y.neg_in_place();
        *y += x;
    }
}