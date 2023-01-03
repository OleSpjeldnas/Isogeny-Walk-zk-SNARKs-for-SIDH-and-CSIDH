def prime_string(n):
  def is_prime(num):
    # 1 is not considered a prime number
    if num == 1:
      return False
    # 2 and 3 are prime numbers
    if num == 2 or num == 3:
      return True
    # Even numbers are not prime, so we can check only odd numbers
    if num % 2 == 0:
      return False
    # Check for prime by checking for divisibility up to the square root of the number
    for i in range(3, int(num ** 0.5) + 1, 2):
      if num % i == 0:
        return False
    # If the number is not divisible by any number up to its square root, it is prime
    return True

  primes = []
  num = 2
  while len(primes) < n:
    if is_prime(num):
      primes.append(num)
    num += 1
  
  return "(" + ", ".join(str(p) for p in primes) + ")"

#print(prime_string(75))

def product_of_first_n_primes(n):

  def is_prime(num):
    # 1 is not considered a prime number
    if num == 1:
      return False
    # 2 and 3 are prime numbers
    if num == 2 or num == 3:
      return True
    # Even numbers are not prime, so we can check only odd numbers
    if num % 2 == 0:
      return False
    # Check for prime by checking for divisibility up to the square root of the number
    for i in range(3, int(num ** 0.5) + 1, 2):
      if num % i == 0:
        return False
    # If the number is not divisible by any number up to its square root, it is prime
    return True

  primes = []
  num = 2
  while len(primes) < n:
    if is_prime(num):
      primes.append(num)
    num += 1
  
  product = 1
  for prime in primes:
    product *= prime
  
  return product

def prime_factors(n):
  def is_prime(num):
    # 1 is not considered a prime number
    if num == 1:
      return False
    # 2 and 3 are prime numbers
    if num == 2 or num == 3:
      return True
    # Even numbers are not prime, so we can check only odd numbers
    if num % 2 == 0:
      return False
    # Check for prime by checking for divisibility up to the square root of the number
    for i in range(3, int(num ** 0.5) + 1, 2):
      if num % i == 0:
        return False
    # If the number is not divisible by any number up to its square root, it is prime
    return True

  factors = []
  # Check for divisibility by 2
  while n % 2 == 0:
    factors.append(2)
    n /= 2
  # Check for divisibility by odd numbers
  for i in range(3, int(n ** 0.5) + 1, 2):
    while n % i == 0:
      factors.append(i)
      n /= i
  # If the number is greater than 2, it must be prime
  if n > 2:
    factors.append(n)
  return factors
def fac_of_prod(n):
    p = product_of_first_n_primes(n) - 1
    return prime_factors(p)
def double_fac_of_prod(n):
    p = 2*product_of_first_n_primes(n) - 1
    return prime_factors(p)

print(1223*(2**213)*(3**133)+1)