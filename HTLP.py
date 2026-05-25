from Crypto.Util.number import getPrime, getStrongPrime, isPrime, getRandomRange, GCD \
, inverse
import gmpy2


def sample_Z_star(N):
    while True:
        s = getRandomRange(1,N)
        if GCD(s,N) == 1:
            return s


class LHP:

    def __init__(self,T,N,g,h):
        self.T = T
        self.N = N
        self.g = g
        self.h = h
        self.pp = (T,N,g,h)

    def PSetup (l,T):
        p = getPrime(l - 1) + 1
        q = getPrime(l - 1) + 1
        N = p*q
        N2 = N*N
        g_tilde = sample_Z_star(N)
        g = (-pow(g_tilde,2,N)) % N
        temp_exp = pow(2,T,(p - 1)*(q - 1)//2) # calcula 2^T modulo phi(N),  para facilitar a proxima operacao
        h = pow(g,temp_exp,N2) # h só vai ser utilizado numa potencia sua mod (N2), fazemos isto para facilitar
        return LHP(T,N,g,h)
    
    def PGen(pp,s):
        N2 = pp.N*pp.N
        r = getRandomRange(1,N2)
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
        inv_wN = inverse(pow(w,N,N2),N2) # facilita a próxima operação
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

pp = LHP.PSetup(512,16)
Z1 = LHP.PGen(pp,12)
Z2 = LHP.PGen(pp,512)
Z_prime = LHP.PEval(pp,(Z1,Z2))
s_prime = LHP.PSolve(pp,Z_prime)
print(s_prime)

        



