from Crypto.Util.number import getPrime, getStrongPrime, isPrime, getRandomRange, GCD
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
        g_tilde = sample_Z_star(N)
        g = (-pow(g_tilde,2,N)) % N
        temp_exp = pow(2,T,(p - 1)*(q - 1)//2) # calcula 2^T modulo phi(N),  para facilitar a proxima operacao
        h = pow(g,temp_exp)
        return LHP(T,N,g,h)
    
    def PGen(pp,s):
        N2 = pp.N*pp.N
        r = getRandomRange(1,N2)
        u = pow(pp.g,r,pp.N)
        v = pow(pp.h,r*pp.N,N2)*pow((1 + pp.N),s,N2)
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

        s = ((v // pow(w,N,N2)) % N2 - 1) // N
        return s                
    

pp = LHP.PSetup(512,21)
Z = LHP.PGen(pp,12)
s = LHP.PSolve(pp,Z)
print(s)

        



