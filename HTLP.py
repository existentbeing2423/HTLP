from Crypto.Util.number import getPrime, getRandomRange, GCD , inverse
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
        p = 2*getPrime(l) + 1
        q = 2*getPrime(l) + 1
        N = p*q
        g_tilde = sample_Z_star(N)
        g = (-pow(g_tilde,2,N)) % N
        temp_exp = pow(2,T,(p - 1)*(q - 1)//2) # calcula 2^T modulo phi(N)/2,  para facilitar a proxima operacao
        h = pow(g,temp_exp,N) 
        print("done setup \n")
        return LHP(T,N,g,h)
    
    def PGen(pp,s):
        N2 = pp.N*pp.N
        r = getRandomRange(1,N2)
        u = pow(pp.g,r,pp.N)
        v = pow(pp.h,r*pp.N,N2)*pow((1 + pp.N),s,N2) % N2
        print("done gen\n")
        return (u,v)
    
    def PSolve(pp,Z):
        T = pp.T
        N = pp.N
        N2 = N*N
        w = Z[0]
        v = Z[1]
        i = 0
        while i < T: # calculo de u^2T por quadratura repetida
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
        print("done eval\n")
        return (u_prime,v_prime)


class MHP:

    def __init__(self,T,N,g,h):
        self.T = T
        self.N = N
        self.g = g
        self.h = h
        self.pp = (T,N,g,h)

    def PSetup (l,T):
        p = 2*getPrime(l) + 1
        q = 2*getPrime(l) + 1
        N = p*q
        g_tilde = sample_Z_star(N)
        g = (-pow(g_tilde,2,N)) % N
        temp_exp = pow(2,T,(p - 1)*(q - 1)//2) # calcula 2^T modulo phi(N),  para facilitar a proxima operacao
        h = pow(g,temp_exp,N) 
        print("done setup \n")
        return LHP(T,N,g,h)
    
    def PGen(pp,s):
        N = pp.N
        N2 = N*N
        r = getRandomRange(1,N2)
        u = pow(pp.g,r,pp.N)
        v = pow(pp.h,r,N)*(s % N) % N
        vh = pow(pp.h,r,N)
        print("done gen\n")
        return (u,v)
    
    def PSolve(pp,Z):
        T = pp.T
        N = pp.N
        w = Z[0]
        v = Z[1]
        i = 0
        while i < T: # calculo de u^2T por quadratura repetida
            w = pow(w,2,N)
            i += 1
        inv_w = inverse(w,N) # facilita a próxima operação
        s = (v * inv_w) % N 
        return s

    def PEval(pp,Z_list):
        N = pp.N
        N2 = N*N
        u_prime = 1
        v_prime = 1
        for Z in Z_list:                
            u_prime = (u_prime % N * Z[0] % N) % N
            v_prime = (v_prime % N * Z[1] % N) % N
        print("done eval\n")
        return (u_prime,v_prime)


pp = MHP.PSetup(1024,2041)
Z1 = MHP.PGen(pp,12)
Z2 = MHP.PGen(pp,512)
Z_prime = MHP.PEval(pp,(Z1,Z2))
s_prime = MHP.PSolve(pp,Z_prime)
print(s_prime)


        



