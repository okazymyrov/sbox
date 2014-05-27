r"""
    Generation functions

AUTHORS:

- Oleksandr Kazymyrov (2011-06-21): initial version

- Oleksandr Kazymyrov (2013-08-04): rewrote lots of functions

"""

#*****************************************************************************
#       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***************************************************************************** 

def gen_APN6(self, **kwargs):
    r"""
    Generate APN function for ``n=m=6``
    
    EXAMPLE::

        sage: S=Sbox(n=6,m=6)
        sage: S.generate_sbox(method='APN6')
        sage: S.is_APN()
        True
        sage: S.MDT()
        2
    
    """
    if self._n != 6:
        raise TypeError("n must be equal 6")

    F="g*x^3+g^5*x^10+g^4*x^24"

    tr=self._P(self.Tr_pol(x=self._P("x"),n=self._n,m=self._n/2))

    #M = matrix(GF(2),2*self._n)
    
    M1 = self.l2m(tr.subs(self._P("(%s)*x"%(self._alpha^4))))
    M2 = self.l2m(self._alpha*tr)
    M3 = self.l2m(tr.subs(self._P("(%s)*x"%(self._alpha))))
    M4 = self.l2m(self._alpha*tr.subs(self._P("(%s)*x"%(self._alpha^4))))
    
    #M[:self._n,:self._n] = M1
    #M[:self._n,self._n:] = M2
    
    #print "M:\n{0}".format(M.echelon_form())

    #M[:self._n,:self._n] = M3
    #M[:self._n,self._n:] = M4
    
    #print "\nM:\n{0}".format(M.echelon_form())
    
    self._CCZ(G=F,M1=M1,M2=M2,M3=M3,M4=M4)

def gen_CCZ(self,**kwargs):
    r"""
        Compute CCZ-equivalence using the form::

            Y = M * X

                [  x   ]     [ M1 | M2 ]    [ F1(x) ]
            X = [------]  M =[---------] Y =[-------]
                [ G(x) ]     [ M3 | M4 ]    [ F2(x) ]

            F = F2(F1_1(x))

        Functions ``F`` and ``G`` are CCZ-equivalent, F1_1 is inverse function to F1.

        If ``G`` is not presented, then ``G`` computes using the substitution.

        INPUT::

            - ``G`` -- input polynomial
            - ``M1`` -- the matrix corresponding to ``L_1``
            - ``M2`` -- the matrix corresponding to ``L_2``
            - ``M3`` -- the matrix corresponding to ``L_3``
            - ``M4`` -- the matrix corresponding to ``L_4``
            - ``M`` -- see the description above

        EXAMPLE::

            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='polynomial',G="x^3",T='CCZ')
            sage: S.is_APN()
            True

            sage: S
            [0, 7, 0, 1, 5, 0, 4, 7]

            sage: S=Sbox(n=6,m=6)
            sage: M1 = matrix(GF(2),6,[[0,1,0,1,1,0],[1,0,0,0,0,0],[0,1,1,1,1,0],[0,0,0,0,0,0],[0,1,1,1,1,0],[1,1,1,1,1,0]])
            sage: M2 = matrix(GF(2),6,[[0,1,0,1,1,1],[0,0,1,0,1,0],[0,0,1,1,1,0],[0,0,1,1,1,0],[0,1,0,1,1,1],[0,1,1,0,0,1]])
            sage: M3 = matrix(GF(2),6,[[1,1,1,0,1,0],[0,1,1,1,0,0],[1,1,0,0,1,1],[0,0,0,0,0,0],[1,1,0,0,1,1],[1,0,1,1,1,1]])
            sage: M4 = matrix(GF(2),6,[[1,1,1,1,1,0],[1,0,1,0,0,0],[1,0,0,0,0,0],[1,0,0,0,0,0],[1,1,1,1,1,0],[0,1,1,1,1,0]])
            sage: S.generate_sbox(method='polynomial',G="g*x^3+g^5*x^10+g^4*x^24",T='CCZ',M1=M1,M2=M2,M3=M3,M4=M4)
            sage: S.is_bijection()
            True

            sage: S.is_APN()
            True
    """

    G = self.g2p(kwargs.get('G',None))
    M = kwargs.get('M','random')

    if G == None:
        if self._S is not None:
            G = self.interpolation_polynomial(representation="internal")
            self._polynomial = None
        elif self._polynomial is not None:
            G = self._polynomial
            self._polynomial = None
        else:
            G = self._P("x")

    tS = [G.subs(self._K.fetch_int(g)).integer_representation() for g in xrange(self._length)]

    if ('M1' in kwargs) or ('M2' in kwargs) or ('M3' in kwargs) or ('M4' in kwargs):
        M1=kwargs.get('M1',identity_matrix(GF(2),self._n))
        M2=kwargs.get('M2',zero_matrix(GF(2),self._m))
        M3=kwargs.get('M3',identity_matrix(GF(2),self._n))
        M4=kwargs.get('M4',zero_matrix(GF(2),self._m))

        M=matrix(GF(2),self._n+self._m,self._n+self._m)

        M.set_block(0,0,M1)
        M.set_block(0,self._n,M2)
        M.set_block(self._n,0,M3)
        M.set_block(self._n,self._n,M4)

    else:
        if M == "random":

            E = identity_matrix(GF(2),self._n)
            Z = zero_matrix(GF(2),self._n)

            while True:
                self._S=[-1 for g in xrange(self._length)]

                M12 = random_matrix(GF(2),self._n,self._n+self._m,algorithm='echelonizable',rank=self._n)

                for x in xrange(self._length):
                    X = vector(ZZ(x).digits(base=2,padto=self._n)+ZZ(tS[x]).digits(base=2,padto=self._m))
                    X = X.column()
                    X = M12*X
                    X = X.list()
                    self._S[x] = ZZ(X[0:self._n],2)

                if (len(set(self._S)) == (self._length)) and (-1 not in self._S):
                    M = zero_matrix(GF(2),2*self._n,self._n+self._m)
                    M.set_block(0,0,M12)

                    while M.is_singular():
                        M34 = random_matrix(GF(2),self._n,self._n+self._m,algorithm='echelonizable',rank=self._n)
                        M.set_block(self._n,0,M34)

                    break

    self._S=[-1 for g in xrange(self._length)]

    for x in xrange(self._length):
        X = vector(ZZ(x).digits(base=2,padto=self._n)+ZZ(tS[x]).digits(base=2,padto=self._m))
        X = X.column()
        X = M*X
        X = X.list()
        self._S[ZZ(X[0:self._n],2)] = ZZ(X[self._n:],2)

    if -1 in self._S:
        self._S = None
        raise TypeError("Sbox has '-1'. Somethink wrong in 'CCZ'!")

    self._polynomial = None

def gen_dickson(self, **kwargs):
    r"""
    Generate a Dickson polynomial, which is a permutation polynomial [LIDL]
    
    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='dickson')
        sage: S # random

        sage: S.is_bijection()
        True

    REFERENCES::

    [LIDL] R. Lidl and H. Niederreiter, Finite Fields, Encyclopedia of Mathematics and its Applications v. 20, pt. 1, Cambridge University Press, 1997. - pp. 355-356
    """

    while (1):
        k=ZZ(randint(2,floor(self._length-2)))
        if gcd(k,self._length*self._length-1) == 1:
            break

    x=self._P.gen()

    self._polynomial=sum(GF(2)((k/(k-j))*binomial(k-j,j))*((-self._alpha)^j)*(x^(k-2*j)) for j in xrange(floor(k/2)+1))

    self._polynomial=self._polynomial.mod(self._P("x^%d+x"%(self._length)))

def gen_dobbertin(self, **kwargs):
    r"""
    Generate Dobbertin's polynomial. ``\delta``-uniformity equals 2 for odd ``n``.
    
    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='dobbertin')
        sage: S.MDT()
        2

    """

    if (self._n % 5 != 0) or (self._m % 5 != 0):
        raise TypeError("m != 0 mod 5")

    h=self._m/5

    self._polynomial=self._P("x^(2^(4*{0})+2^(3*{0})+2^(2*{0})+2^({0})-1)".format(h))

def gen_EA(self, **kwargs):
    r"""
    Compute ``F(x) = M1 * G( M2 * x + V2 ) + M3 * x + V1``.

    If ``G`` is not presented, then ``G`` computes using the substitution.

    INPUT::

        - ``G`` -- input polynomial
        - ``M1, M2`` -- two linear permutation polynomials given in matrix form
        - ``M3`` -- arbitrary linear polynomials
        - ``V1`` -- two vectors in ``F_2^m``
        - ``V2`` -- two vectors in ``F_2^n``

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^3',T='EA')
        sage: S.MDT()
        2

        sage: S.get_linear_functions() # random
        [
        [1 0 1 0 0]  [0 1 0 1 1]  [0 1 1 1 1]                 
        [0 1 1 0 0]  [1 0 0 1 0]  [1 0 0 1 0]                 
        [0 0 0 0 1]  [1 0 1 1 0]  [1 1 1 0 1]                 
        [0 0 1 1 1]  [0 0 0 0 1]  [1 1 0 1 0]                 
        [0 0 0 1 1], [0 0 1 1 0], [0 1 0 0 1], (1, 0, 0, 1, 0),

        (0, 0, 0, 1, 0)
        ]

        sage: S.generate_sbox(method='polynomial',G='x^3',T='A')
        sage: S.MDT()
        2

        sage: S.is_bijection()
        True

        sage: S.generate_sbox(method='polynomial',G='x^3',T='EA',M1=random_matrix(GF(2),5,5,algorithm='echelonizable',rank=5))
        sage: S.get_linear_functions() # random
        [
        [0 0 1 1 1]                        
        [1 0 0 1 0]                        
        [1 1 1 1 0]                        
        [1 1 1 0 0]                        
        [0 0 0 0 1], None, None, None, None
        ]

        sage: S.is_APN()
        True

        sage: S.is_bijection()
        True
    """
    G = self.g2p(kwargs.get('G',None))
    
    M1 = kwargs.get('M1',None)
    V1 = kwargs.get('V1',None)
    M2 = kwargs.get('M2',None)
    V2 = kwargs.get('V2',None)
    M3 = kwargs.get('M3',None)
    
    if G == None:
        if self._S is not None:
            G = self.interpolation_polynomial(representation="internal")
            self._polynomial = None
        elif self._polynomial is not None:
            G = self._polynomial
            self._polynomial = None
        else:
            G = self._P("x")

    self._S=range(self._length)

    if M1 == "random":
        M1=random_matrix(GF(2),self._m,self._m)
        while(M1.is_singular()):
            M1=random_matrix(GF(2),self._m,self._m)
    elif M1 is not None:
        M1=Matrix(GF(2),M1)

    if M2 == "random":
        M2=random_matrix(GF(2),self._n,self._n)
        while(M2.is_singular()):
            M2=random_matrix(GF(2),self._n,self._n)
    elif M2 is not None:
        M2=Matrix(GF(2),M2)

    if M3 == "random":
        M3=random_matrix(GF(2),self._n,self._m)
    elif M3 is not None:
        M3=Matrix(GF(2),M3)

    if V1 == "random":
        V1=random_vector(GF(2),self._m)
    elif V1 is not None:
        V1=vector(GF(2),V1)    

    if V2 == "random":
        V2=random_vector(GF(2),self._n)
    elif V2 is not None:
        V2=vector(GF(2),V2)

    #print "~"*40
    #print "M1:\n{0}".format(M1)
    #print "V1:\n{0}".format(V1)
    #print "M2:\n{0}".format(M2)
    #print "V2:\n{0}".format(V2)
    #print "M3:\n{0}".format(M3)
    #print "~"*40

    self._LFoEA = [M1,M2,M3,V1,V2]

    for i in xrange(self._length):
        self._S[i]=vector(GF(2),ZZ(i).digits(base=2,padto=self._n))

        if M2 is not None:
            self._S[i]=M2*self._S[i]

        if V2 is not None:
            self._S[i]=vector(GF(2),[ZZ(self._S[i][j]) ^^ ZZ(V2.get(j)) for j in xrange(len(self._S[i]))])

        self._S[i]=self._K(PolynomialRing(GF(2),'x')(self._S[i].list()).mod(self._K.modulus()).list())
        self._S[i]=G.subs(self._S[i])
        self._S[i]=self._S[i]._vector_()

        if M1 is not None:
            self._S[i]=M1*self._S[i]

        if M3 is not None:
            tx=M3*vector(GF(2),ZZ(i).digits(base=2,padto=self._m))
            self._S[i]=vector(GF(2),[ZZ(self._S[i].get(j)) ^^ ZZ(tx.get(j)) for j in xrange(len(self._S[i]))])

        if V1 is not None:
            self._S[i]=vector(GF(2),[ZZ(self._S[i].get(j)) ^^ ZZ(V1.get(j)) for j in xrange(len(self._S[i]))])

        self._S[i]=ZZ(self._S[i].list(),2)

def gen_gold(self, **kwargs):
    r"""
    Generate Gold's polynomial. ``\delta``-uniformity equals 2 for odd ``n``.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='gold')
        sage: S.is_APN()
        True
    """

    while (1):
        i=randint(1,floor(self._m/2))
        if gcd(i,self._m) == 1:
            break

    self._polynomial=self._P("x^(2^{0}+1)".format(i))

def gen_inverse(self, **kwargs):
    r"""
    Generate inverse polynomial. ``\delta``-uniformity equals 2 for odd ``n``.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='inverse')
        sage: S.is_APN()
        True
    """
    self._polynomial=self._P("x^(2^({0}-1)-1)".format(self._m))
    
def gen_kasami(self, **kwargs):
    r"""
    Generate Kasami's polynomial. ``\delta``-uniformity equals 2 for odd ``n``.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='kasami')
        sage: S.is_APN()
        True
    """

    while (1):
        i=randint(1,floor(self._m/2))
        if gcd(i,self._m) == 1:
            break

    self._polynomial=self._P("x^(2^(2*{0})-2^{0}+1)".format(i))

def gen_niho(self, **kwargs):
    r"""
    Generate Niho's polynomial. ``\delta``-uniformity equals 2 for odd ``n``.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='niho')
        sage: S.is_APN()
        True

    """

    if is_odd(self._m):
        h=(self._m-1)/2
    else:
        h=(self._m)/2

    if is_even(h):
        self._polynomial=self._P("x^(2^({0})+2^({0}/2)-1)".format(h))
    else:
        self._polynomial=self._P("x^(2^({0})+2^((3*{0}+1)/2)-1)".format(h))

def gen_OP4(self, **kwargs):
    r"""
    Generate optimal permutation polynomial (with minimum degree 3) for ``n=m=4`` [LePo07]

    INPUT::

        - ``function`` --  if presented then return a fixed function (input value from 0 to 7)

    EXAMPLE::

        sage: S=Sbox(n=4,m=4)
        sage: S.generate_sbox(method='OP4')
        sage: S.MDT()
        4

        sage: S.minimum_degree()
        3

        sage: S.generate_sbox(method='OP4',function=1)
        sage: F=S.get_polynomial()
        sage: F
        'x^14 + g^11*x^13 + g^7*x^12 + g*x^11 + g^8*x^10 + g^13*x^9 + g^11*x^8 + g^2*x^6 + g*x^5 + g^2*x^4 + g^7*x^3 + g*x^2 + g^8*x'

        sage: S.generate_sbox(method='OP4',function=1)
        sage: F == S.interpolation_polynomial()
        True

    REFERENCES::

    .. [LePo07] G. Leander and A. Poschmann. "On the Classification of 4 Bit S-boxes." Arithmetic of Finite Fields. Springer Berlin Heidelberg, 2007. 159-176.

    """

    if self._n != 4:
        raise TypeError("Implemented only for n=4")

    G = [
        self.g2p("x^14+g^11*x^13+g*x^12+g^3*x^11+g^5*x^9+g^7*x^8+g^8*x^7+g^4*x^6+g^11*x^5+g^2*x^4+g^4*x^3+g^11*x^2"),
        self.g2p("x^14+g^11*x^13+g^7*x^12+g*x^11+g^8*x^10+g^13*x^9+g^11*x^8+g^2*x^6+g*x^5+g^2*x^4+g^7*x^3+g*x^2+g^8*x"),
        self.g2p("x^14+g^13*x^13+g^9*x^12+g^6*x^11+g^10*x^10+g^7*x^9+g^10*x^8+g^7*x^7+g^8*x^6+g^12*x^5+g^12*x^4+x^3+g^11*x^2+g^5*x"),
        self.g2p("x^14+g^4*x^13+g^3*x^12+g^2*x^11+x^10+g^11*x^9+g^2*x^8+g*x^7+g^2*x^6+g^9*x^5+g^4*x^4+g^9*x^3+g^12*x^2+g^11*x"),
        self.g2p("x^14+g*x^13+g^9*x^12+g*x^11+g^7*x^10+g^6*x^7+g^10*x^6+g*x^5+g^8*x^4+g^2*x^3+g^6*x^2+g^9*x"),
        self.g2p("x^14+g*x^13+x^12+g^7*x^11+g^13*x^10+g*x^9+g^11*x^8+g^14*x^7+g^3*x^6+g^6*x^5+g*x^4+g^14*x^3+g^14*x^2+g^9*x"),
        self.g2p("x^14+g^10*x^13+g*x^12+g^4*x^11+g^14*x^10+g^4*x^9+g^5*x^8+g^2*x^7+g^9*x^6+g^4*x^5+g^8*x^4+g^14*x^3+g^5*x^2+x"),
        self.g2p("x^14+g^12*x^13+g^8*x^12+g^8*x^11+g^14*x^10+g*x^9+g^8*x^8+g^14*x^7+g^6*x^6+x^5+g^14*x^4+g^12*x^3+g*x^2+g^14*x")
    ]

    self._polynomial = G[kwargs.get('function',randint(0,7))]

def gen_random_substitution(self, **kwargs):
    r"""
    Generate a random substitution

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='random_substitution')
        sage: S # random
        [2, 4, 7, 7, 4, 5, 2, 7]

    """
    self._S = [randint(0,(1<<self._m)-1) for _ in xrange(self._length)]

def gen_random_permutation(self, **kwargs):
    r"""
    Generate a random permutation

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='random_permutation')
        sage: S.is_bijection()
        True
    """

    if self._n == self._m:
        self._S = range(self._length)
        shuffle(self._S)
    else:
        raise TypeError("m doesn't equal n")

def gen_welch(self, **kwargs):
    r"""
    Generate Welch's polynomial. ``\delta``-uniformity equals 2 for odd ``n``.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='welch')
        sage: S.is_APN()
        True

    """

    if is_odd(self._m):
        h=(self._m-1)/2
    else:
        h=(self._m)/2

    self._polynomial=self._P("x^(2^{0}+3)".format(h))
