import secrets
import time

def gcd(a: int, b: int) -> int:
    while b != 0:
        a, b = b, a % b
    return abs(a)

def sample_Z_star(N):
    while True:
        s = secrets.randbelow(N-1) + 1  # returns a random integer in the range [1, N-1]
        if gcd(s,N) == 1:
            return s


def generate_prime(bits: int) -> int:
    """
    Generate a prime number with the specified number of bits.
    Is forcing to be odd by setting the least significant bit to 1, and ensuring the most significant bit is 1 to get the desired bit length.
    """
    while True:
        candidate = secrets.randbits(bits) | (1 << (bits - 1)) | 1  # OR's forces the first bit to 1 to get higher prime and last bit to be 1 to be odd 
        if miller_rabin(candidate):
            return candidate

def miller_rabin(n: int, k: int = 40) -> bool:
    """
    Test if a number is prime using the Miller-Rabin probabilistic test.
    k is the number of rounds (40 rounds gives a negligible error probability).
    """
    if n == 2 or n == 3:
        return True
    if n < 2 or n % 2 == 0:
        return False

    # Write n-1 as 2^r * d
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Witness loop
    for _ in range(k):
    
        a = secrets.randbelow(n - 3) + 2    #Note that secrets.randbelow(n - 3) returns a random integer in the range [0, n-4], so adding 2 shifts it to [2, n-2].
        
        x = pow(a, d, n) 
        
        if x == 1 or x == n - 1:
            continue
            
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False  # Definitely not prime
            
    return True  # Probably prime

class LHP:

    def __init__(self,T,N,g,h):
        self.T = T
        self.N = N
        self.g = g
        self.h = h
        self.pp = (T,N,g,h)

    def PSetup (l,T):
        
        p = generate_prime(l) 
        print(p)
        q = generate_prime(l)
        print(q)
        
        N = p*q
        N2 = N*N
        g_tilde = sample_Z_star(N)
        g = (-pow(g_tilde,2,N)) % N
        temp_exp = pow(2,T,(p - 1)*(q - 1)//2) # calcula 2^T modulo phi(N),  para facilitar a proxima operacao
        h = pow(g,temp_exp,N) # h só vai ser utilizado numa potencia sua mod (N2), fazemos isto para facilitar
        return LHP(T,N,g,h)
    
    def PGen(pp,s):
        N2 = pp.N*pp.N
        r = secrets.randbelow(N2) + 1 # secrets.randbelow(N2) returns 0 to N^2 - 1. Adding 1 shifts it to 1 to N^2.
        u = pow(pp.g,r,pp.N)
        v = pow(pp.h,r*pp.N,N2)*pow((1 + pp.N),s,N2) % N2
        return (u,v)
    
    def PSolve(pp,Z):
        T = pp.T
        N = pp.N
        N2 = N*N
        w = Z[0]
        v = Z[1]
        i = 0
        while i < T:
            w = pow(w,2,N)
            i += 1
        wN = pow(w, N, N2)
        inv_wN = pow(wN, -1, N2) # facilita a próxima operação
        s = ((v * inv_wN) % N2 - 1) // N
        return s

    def PEval(pp,Z_list):
        N = pp.N
        N2 = N*N
        u_prime = 1
        v_prime = 1
        for Z in Z_list:                
            u_prime = (u_prime % N * Z[0] % N) % N
            v_prime = (v_prime % N2 * Z[1] % N2) % N2
        return (u_prime,v_prime)

t_start = time.time()

pp = LHP.PSetup(2048,16)
Z1 = LHP.PGen(pp,142)
Z2 = LHP.PGen(pp,512)
Z_prime = LHP.PEval(pp,(Z1,Z2))
s_prime = LHP.PSolve(pp,Z_prime)
print(s_prime)

t_end = time.time()
print(f"  [Solve] Done in {t_end - t_start:.4f}s")