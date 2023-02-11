use ark_ff::{Fp512, MontBackend, Fp2, MontFp, Fp2Config};


#[derive(ark_ff::fp::MontConfig)]
#[modulus = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659"]
#[generator = "2"]

pub struct FqConfig;
pub type Kp = Fp512<MontBackend<FqConfig, 8>>;

pub type K = Fp2<Fq2Config>;
pub struct Fq2Config;
impl Fp2Config for Fq2Config {
    type Fp = Kp;

    /// NONRESIDUE = -1
    const NONRESIDUE: Kp = MontFp!("-1");

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP2_C1: &'static [Kp] = &[
        MontFp!("-1"),
        MontFp!("-1")
    ];

}