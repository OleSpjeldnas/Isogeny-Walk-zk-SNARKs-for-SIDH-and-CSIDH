use ark_ff::{Fp512, MontBackend, BitIteratorBE, One};
use core::ops::{MulAssign,Mul};
#[derive(ark_ff::fp::MontConfig)]
#[modulus = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659"]
#[generator = "2"]

pub struct FqConfig;
pub type F = Fp512<MontBackend<FqConfig, 8>>;

trait Pow:
'static
+ One
+ Mul<Self, Output = Self>
+ MulAssign<Self>
+ for<'a> Mul<&'a Self, Output = Self>
+ for<'a> MulAssign<&'a mut Self>
+ for<'a> MulAssign<&'a Self> 
{
    fn square_in_place(&mut self) -> &mut Self;

    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one();

        for i in BitIteratorBE::without_leading_zeros(exp) {
            res.square_in_place();

            if i {
                res *= self;
            }
        }
        res
    }
}

impl Pow for F {
    fn square_in_place(&mut self) -> &mut Self {
        self.square_in_place()
    }
}

