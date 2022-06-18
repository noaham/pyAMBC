from ambc import *
from itertools import product

class Perm(object):
    def __init__(self, wind):
        self.window = wind
        self.n = len(self.window)

    def __call__(self,i):
        # returns the effect of the permutation on i
        # i = r + kn for some r = 1,2,...,n
        r = modn(i,self.n)
        return self.window[r-1]+i-r

    def __eq__(self,other):
        return self.window == other.window
        
    def __mul__(self,other):
        if self.n == other.n:
            return Perm([self(other(i)) for i in range(1, self.n + 1)])
        else:
            raise ValueError("permutations must be the same rank")

    def __repr__(self):
        return str(self.window)

    def inv(self):
        winv = ['x']*self.n
        for i in range(1,self.n+1):
            wi = self(i)
            # say wi = r + kn, we want r = wi - kn --> i - kn 
            r = modn(wi,self.n)
            winv[r-1] = i - (wi-r)
        return Perm(winv)

    def __len__(self):
        # returns the length of the permutation
        len = 0
        for i in range(1,self.n+1):
            for j in range(i+1,self.n+i+1):
                t = 0
                while self(j+t*self.n) < self(i):
                    t += 1
                    len += 1
        return len
        
    def ldesc(self):
        # returns the left descent set of
        pinv = self.inv()
        return [i for i in range(self.n) if pinv(i) > pinv(i+1)]

    def rdesc(self):
        # returns the left descent set of
        return [i for i in range(self.n) if self(i) > self(i+1)]

    # the following two functions define permutations w = ax such that x is
    # in the parabolic subgroup permuting 1,2,...,n and a is the minimal
    # coset representative
    def a(self):
        return Perm(sorted(self.window))
    
    def x(self):
        return self.a().inv()*self

    def AMBC(self):
        return AMBC(self.window)

    def RSK(self):
        if self.a() == id(self.n):
            return RSK(self.window)
        else:
            raise ValueError("permutation must be in the finite symmetric group")

    def cactus(self):
        x = self.x()
        a = self.a()
        (P,Q) = x.RSK()
        eQ = Q.schutzenberger()
        ex = Perm(invRSK(P, eQ))
        ew = a*ex
        return ew

    def schutz(self):
        [P,Q,rho] = self.AMBC()
        return Perm(invAMBC(P, Q.schutzenberger(), rho))

    def rev_comp(self):
        return Perm([self.n - w + 1 for w in self.window[::-1]])
            
def s(i,n):
    # the permutation that swaps i, i+1
    if i == 0:
        return Perm([0] + [i for i in range(2, n)] + [n + 1])
    elif i in range(1,n):
        return Perm([i for i in range(1, i)] + [i + 1, i] + [i for i in range(i + 2, n + 1)])
    else:
        raise ValueError("i must be between 0 and n-1")

def id(n):
    return Perm(list(range(1, n + 1)))

def prod(word,n=None):
    if word == []:
        if n == None:
            raise ValueError("if word is empty, need to specify n")
        return id(n)
    product = id(word[0].n)
    for p in word:
        product = product*p
    return product

import random
def rand_perms(N,n=None,type='affine'):
    # returns a list of N random affine permutations
    # if n=None, then they are random amongst the groups hatS_3 to hatS_9
    # a random perm is chosen by randomly choosing a integer l (from 1 to 100)
    # and randomly choosing l simple reflections and multiplying them together.
    perms = []
    for i in range(N):
        if n == None:
            m = random.randint(3,9)
        else:
            m = n
        l = random.randint(1,100)
        if type == 'affine':
            t=0
        elif type == 'finite':
            t=1
        word = [s(random.randint(t,m-1),m) for i in range(l)]
        perms.append(prod(word))
    return perms

def conj_test(fun,N,n=None,type='affine'):
    perms = rand_perms(N,n,type)
    return all([fun(p) for p in perms])

def conj_testr(fun,N,n=None,type='affine'):
    perms = rand_perms(N,n,type)
    for p in perms:
        if not fun(p):
            return p
    return all([fun(p) for p in perms])
