from Crypto.Util.number import getRandomRange, GCD , inverse
import time 
from matplotlib import pyplot as plt
import subprocess as sp


def sample_Z_star(N):
    while True:
        s = getRandomRange(1,N)
        if GCD(s,N) == 1:
            return s


class Partial_HP:

    def __init__(self,T,N,g,h):
        self.T = T
        self.N = N
        self.g = g
        self.h = h
        self.pp = (T,N,g,h)

    def PSetup (l,T):
        while True:
            p = int(
                sp.check_output(
                    ["openssl", "prime", "-generate", "-safe", "-bits", "1024"]
                )
            )
            q = int(
                sp.check_output(
                    ["openssl", "prime", "-generate", "-safe", "-bits", "1024"]
                )
            )
            if p != q:
                break
        N = p*q
        g_tilde = sample_Z_star(N)
        g = (-pow(g_tilde,2,N)) % N
        temp_exp = pow(2,T,(p - 1)*(q - 1)//2) # calcula 2^T modulo phi(N)/2,  para facilitar a proxima operacao
        h = pow(g,temp_exp,N) 
        return Partial_HP(T,N,g,h)
    
    def LPGen(pp,s):
        N2 = pp.N*pp.N
        r = getRandomRange(1,N2)
        u = pow(pp.g,r,pp.N)
        v = pow(pp.h,r*pp.N,N2)*pow((1 + pp.N),s,N2) % N2
        return (u,v)
    
    def LPSolve(pp,Z):
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

    def LPEval(pp,Z_list):
        N = pp.N
        N2 = N*N
        u_prime = 1
        v_prime = 1
        for Z in Z_list:                
            u_prime = (u_prime % N * Z[0] % N) % N
            v_prime = (v_prime % N2 * Z[1] % N2) % N2
        return (u_prime,v_prime)
    
    def MPGen(pp,s):
        N = pp.N
        N2 = N*N
        r = getRandomRange(1,N2)
        u = pow(pp.g,r,pp.N)
        v = pow(pp.h,r,N)*(s % N) % N
        return (u,v)
    
    def MPSolve(pp,Z):
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

    def MPEval(pp,Z_list):
        N = pp.N
        u_prime = 1
        v_prime = 1
        for Z in Z_list:                
            u_prime = (u_prime % N * Z[0] % N) % N
            v_prime = (v_prime % N * Z[1] % N) % N
        return (u_prime,v_prime)



T = 100
T_axis = []
tsetup_axis = []
tgen_axis = []
teval_axis = []
tsolve_axis = []

while T < 10000000:

    T_axis.append(T)

    tsetup_start = time.time()
    pp = Partial_HP.PSetup(1024,T)
    tsetup_end = time.time()
    tsetup_axis.append(tsetup_end - tsetup_start)

    tgen_start = time.time()
    Z1 = Partial_HP.MPGen(pp,12)
    Z2 = Partial_HP.MPGen(pp,512)
    tgen_end = time.time()
    tgen_axis.append(tgen_end - tgen_start)


    teval_start = time.time()
    Z_prime = Partial_HP.MPEval(pp,(Z1,Z2))
    teval_end = time.time()
    teval_axis.append(teval_end - teval_start)


    tsolve_start = time.time()
    s_prime = Partial_HP.MPSolve(pp,Z_prime)
    tsolve_end = time.time()
    tsolve_axis.append(tsolve_end - tsolve_start)

    if s_prime != 512 * 12:
        raise Exception("Error: Wrong result")

    T = T * 2

plt.figure()
plt.plot(T_axis,tsetup_axis,label = "Setup time")
plt.plot(T_axis,tgen_axis,label = "Gen time")
plt.plot(T_axis,tsolve_axis,label = "Solve time")
plt.xscale('log')
plt.xlabel("T (time complexity factor)")
plt.ylabel("time")
plt.title("Multiplicatively HTLP")
plt.legend()
plt.savefig('MHTLP.png')
print("s' = s1 * s2 = ",s_prime)

T = 100
T_axis = []
tsetup_axis = []
tgen_axis = []
teval_axis = []
tsolve_axis = []

while T < 10000000:

    T_axis.append(T)

    tsetup_start = time.time()
    pp = Partial_HP.PSetup(1024,T)
    tsetup_end = time.time()
    tsetup_axis.append(tsetup_end - tsetup_start)


    tgen_start = time.time()
    Z1 = Partial_HP.LPGen(pp,12)
    Z2 = Partial_HP.LPGen(pp,512)
    tgen_end = time.time()
    tgen_axis.append(tgen_end - tgen_start)


    teval_start = time.time()
    Z_prime = Partial_HP.LPEval(pp,(Z1,Z2))
    teval_end = time.time()
    teval_axis.append(teval_end - teval_start)


    tsolve_start = time.time()
    s_prime = Partial_HP.LPSolve(pp,Z_prime)
    tsolve_end = time.time()
    tsolve_axis.append(tsolve_end - tsolve_start)


    if s_prime != 512 + 12:
        raise Exception("Error: Wrong result")

    T = T * 2

plt.figure()
plt.plot(T_axis,tsetup_axis,label = "Setup time")
plt.plot(T_axis,tgen_axis,label = "Gen time")
plt.plot(T_axis,tsolve_axis,label = "Solve time")
plt.xscale('log')
plt.xlabel("T (time complexity factor)")
plt.ylabel("time")
plt.title("Linearly HTLP")
plt.legend()
plt.savefig('LHTLP.png')
print("s' = s1 * s2 = ",s_prime)


        



