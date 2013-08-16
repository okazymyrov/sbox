r"""
    Generation and cryptanalysis of substitutions

    Generation and cryptanalysis function defined externally
    in ``GSbox.sage`` and ``CSbox.sage`` respectively.

AUTHORS::

- Oleksandr Kazymyrov (2011): initial version

- Maksim Storetvedt (2013-04-26): rewrote the function "cycles" in C

- Oleksandr Kazymyrov (2013-08-04): almost complete rewrite

- Anna Maria Eilertsen, Oleksandr Kazymyrov (2013-08-14): added documentations

EXAMPLES::

    sage: S=Sbox(n=3,m=3)
    sage: S.generate_sbox(method="random_permutation")
    sage: S # random
    [1, 2, 3, 5, 0, 4, 7, 6]

"""

#******************************************************************************
#       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
#                     2013 Maksim Storetvedt <Maksim.Storetvedt@student.uib.no>
#                     2013 Anna Maria Eilertsen <Anna.Eilertsen@student.uib.no>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

os.environ['SBOX_ROOT'] = os.getcwd()+"/.."
os.chdir(os.getcwd()+"/..")

load os.environ['SBOX_ROOT'] + "/Cython/CFunc.spyx"
load os.environ['SBOX_ROOT'] + "/Cython/CPPFunc.spyx"
load os.environ['SBOX_ROOT'] + "/Sage/GSbox.sage"
load os.environ['SBOX_ROOT'] + "/Sage/CSbox.sage"

class Sbox(SageObject):
    r"""
    A substitution box or S-box is one of the basic components of
    symmetric key cryptography. In general, an S-box takes ``n`` input
    bits and transforms them into ``m`` output bits. This is called an
    ``nxm`` S-box and is often implemented as a lookup table. 
    
    Another representation is by vectorial Boolean functions, which is defined
    by a Polynomial Ring in ``x`` over Finite Field in ``a`` of size ``2^n`` with 
    primitive element ``g``.
    
    EXAMPLE:
    
    We consider the S-box of the block cipher AES [AES]_::

        sage: sbox=[1,2,3,0,7,6,5,4]
        sage: S=Sbox(n=3,m=3,sbox=[1,2,3,0,7,6,5,4])
        sage: S.get_sbox()
        [1, 2, 3, 0, 7, 6, 5, 4]

    
    Note that by default bits are interpreted in little endian
    notation.
    
    REFERENCES:

        .. [AES] FIPS, PUB. "197: Specification for the Advanced Encryption
        Standard, 2001." http://www.csrc.nist.gov/publications/fips/fips197/fips-197.pdf
        4 (2009): 17-18.

    """

    def __init__(self, n=None, m=None, **kwargs):
        r"""
        
        Construct a substitution box (S-box) with ``n`` input bits and ``m`` output bits 
        with given properties.
        
        If ``n`` and ``m`` are not defined then a ``TypeError`` is raised.

        INPUT::
        
            - ``n``    -- number of input bits of Sbox
        
            - ``m`` -- number of output bits of Sbox
            
            - ``modulus`` -- irreducible polynomial, defining the finite field

            - ``sbox`` -- array of integers, defining a mapping
        
        EXAMPLE::
        
        We construct an ``3x3`` S-box :
::        
            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='random_permutation')
            sage: S # random
            [0, 1, 2, 5, 7, 3, 4, 6]
        
        We construct an ``3x3`` S-box defined by polynomial :
::
            sage: S=Sbox(n=3,m=3,polynomial='g*x^3+1')
            sage: S
            [1, 3, 7, 2, 0, 6, 4, 5]

        We construct an ``3x3`` S-box defined a lookup table :
::
            sage: S=Sbox(n=3,m=3,sbox=[0,2,3,4,1,1,1,1])
            sage: S
            [0, 2, 3, 4, 1, 1, 1, 1]

        Check the correctness of the given arguments  

::
            sage: S=Sbox(n=2,m=2,sbox=[3,9,0,0])
            Traceback (most recent call last):
            ...
            TypeError: The bit length (4) of the maximum element of given S-box larger then internal 'm' (2)

            sage: S=Sbox(n=2,m=2,sbox=[1,2,0,0,0,6])
            Traceback (most recent call last):
            ...
            TypeError: The length (6) of given S-box differs from internal (4)
        """

        if n is None:
            raise TypeError("n is not defined")

        if m is None:
            raise TypeError("m is not defined")

        self._modulus = kwargs.get('modulus',None)
        self._n=n
        self._m=m
        self._length=1<<self._n
        self._LFoEA = [None,None,None,None,None]
        
        if self._modulus is None:
            self._K = GF(1<<max(self._m,self._n),'a',modulus='conway')
        else:
            self._K = GF(1<<max(self._m,self._n),'a',modulus=ZZ['x'](self._modulus))

        self._alpha = self._K.multiplicative_generator()

        self._P = PolynomialRing(self._K,'x')

        self._S = kwargs.get('sbox',None)

        if self._S is not None:
            if len(self._S) != self._length:
                raise TypeError("The length ({0}) of given S-box differs from internal ({1})".format(len(self._S),self._length))

            if ZZ(max(self._S)).nbits() > self._m:
                raise TypeError("The bit length ({0}) of the maximum element of given S-box larger then internal 'm' ({1})".format(ZZ(max(self._S)).nbits(),self._m))

            self._polynomial = None
        else:
            self._polynomial = kwargs.get('polynomial',None)

            if self._polynomial is not None:
                self._EA(G=self._polynomial)

    def __call__(self, val):
        r"""
            Apply substitution to ``val``
            
        INPUT::
            
            -``val`` -- an integer
            
        EXAMPLE::

            sage: S = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: S[1]
            5

            sage: S(4)
            0
        """
        return self._S[val]

    def __cmp__(self, other):
        """
        S-boxes are considered to be equal if all construction
        parameters match.

        EXAMPLE::

            sage: S1 = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: S2 = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: S1 == S2
            True

            sage: S2 = Sbox(n=3,m=3,sbox=[0,0,6,7,0,1,2,3])
            sage: S1 == S2
            False

            # sage: loads(dumps(S1)) == S1 # doesn't work separately
            # True
        """
        return cmp(self._S, other._S)

    def __getitem__(self, val):
        """
        See  :meth:`SBox.__call__`.

        EXAMPLE::

            sage: S = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: S[7]
            3

            sage: S(0)
            4
        """
        return self(val)

    def __iter__(self):
        """
        EXAMPLE::

            sage: S = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: [i for i in S]
            [4, 5, 6, 7, 0, 1, 2, 3]
        """
        for i in xrange(self._length):
            yield self(i)

    def _repr_(self):
        r"""
        EXAMPLE::

            sage: S=Sbox(n=4,m=3)
            sage: S
            []

            sage: S.generate_sbox(method='random_substitution')
            sage: S
            [0, 4, 0, 2, 4, 1, 5, 0, 7, 4, 7, 6, 1, 6, 6, 3]
        """
        if self._S is None:
            return "[]"
        return "[" + ", ".join(map(str,self._S)) + "]"

    def g2p(self,pol=None):
        r"""
        Convert polynomial with primitive element to internal representation
            
        INPUT::

            - ``pol`` -- polynomial
        
        OUTPUT::

            - The internal representation of a polynomial with primitive element
        
        EXAMPLE::

            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='random_permutation')
            sage: S.interpolation_polynomial()
            'g^2*x^5 + g^6*x^3 + g^3*x + g^4'

            sage: S.g2p(S.interpolation_polynomial())
            a^2*x^5 + (a^2 + 1)*x^3 + (a + 1)*x + a^2 + a
        """    
        if pol is not None:
            return self._P(pol.replace("g","({0})".format(self._K.multiplicative_generator())))
        else:
            return None

    def generate_sbox(self, method=None, T=None, **kwargs):
        r'''
        General function for generating substitutions
        
        INPUT::
        
            - ``method`` -- a string, which method will be used, e.g., 
            
                - ``'APN6'`` -- generate APN-permutation function for dimension 6 (m=n=6)
                - ``'dobbertin'`` -- generate substitution based on Dobbertin's polynomial
                - ``'dickson'`` -- generate substitution based on Dobbertin's polynomial
                - ``'gold'`` --  generate substitution based on Gold's polynomial
                - ``'inverse'`` -- generate substitution based on inverse polynomial
                - ``'kasami'`` -- generate substitution based on Kasami's polynomial
                - ``'niho'`` -- generate substitution based on Niho's polynomial
                - ``'polynomial'`` -- generate substitution based on input polynomial
                - ``'random_permutation'`` -- generate random permutation
                - ``'random_substitution'`` -- generate random substitution
                - ``'welch'`` -- generate substitution based on Welch's polynomial

            - ``T`` -- a string, describe an equivalence transformation to apply, e.g. 
                - ``'A'``    -- affine equivalence
                - ``'EA'`` -- extended affine equivalence
                - ``'CCZ'`` -- Carlet-Charpin-Zinoviev equivalence
                - ``'L'``    -- linear equivalence
                
            - ``G`` -- a polynomial for 'polynomial' method
        
        EXAMPLE:: 
        
            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='polynomial',G="g*x^2+x",T='A')
            sage: S # random
            [7, 4, 1, 2, 7, 4, 1, 2]

            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='welch',T='EA')
            sage: S # random
            [3, 7, 6, 1, 0, 1, 1, 3]

            sage: S=Sbox(n=3,m=1)
            sage: S.generate_sbox(method='polynomial',G="x^3")
            sage: S
            [0, 1, 3, 4, 5, 6, 7, 2]

            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='inverse')
            sage: S
            [0, 1, 3, 4, 5, 6, 7, 2]
        '''
        self._polynomial = None
        self._S = None
        self._LFoEA = [None,None,None,None,None]

        methods=dict(AI=self._AI, # experemental
                APN6=self._APN6,
                claude_matrix=self._claude_matrix, # experemental
                dobbertin=self._dobbertin,
                dickson=self._dickson,
                gold=self._gold,
                inverse=self._inverse,
                kasami=self._kasami,
                lilia=self._lilia, # experemental
                lilia_base=self._lilia_base, # experemental
                niho=self._niho,
                OP4=self._OP4,
                polynomial=self._EA,
                random_permutation=self._random_permutation,
                random_substitution=self._random_substitution,
                welch=self._welch
                )

        if not method in methods:
            raise TypeError("Unsupported method '{0}'".format(method))

        if method != 'polynomial':
            methods[method](**kwargs)

        if method == 'random_permutation' or method == 'random_substitution':
            return

        if T == 'CCZ':
            self._CCZ(**kwargs)
            return
        elif T == 'EA':
            if kwargs.has_key('M1') is not True:
                if kwargs.has_key('M2') is not True:
                    if kwargs.has_key('M3') is not True:
                        if kwargs.has_key('V1') is not True:
                            if kwargs.has_key('V2') is not True:
                                kwargs['M1'] = 'random'
                                kwargs['M2'] = 'random'
                                kwargs['M3'] = 'random'
                                kwargs['V1'] = 'random'
                                kwargs['V2'] = 'random'
        elif T == 'A':
            if kwargs.has_key('M1') is not True:
                if kwargs.has_key('M2') is not True:
                    if kwargs.has_key('V1') is not True:
                        if kwargs.has_key('V2') is not True:
                            kwargs['M1'] = 'random'
                            kwargs['M2'] = 'random'
                            kwargs['V1'] = 'random'
                            kwargs['V2'] = 'random'
        elif T == 'L':
            if kwargs.has_key('M1') is not True:
                if kwargs.has_key('M2') is not True:
                    kwargs['M1'] = 'random'
                    kwargs['M2'] = 'random'

        self._EA(**kwargs)
        
    def get_item(self, x):
        r"""
        See `Sbox.__call__`.
        
        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: S(7) == S.get_item(7)
            True
        """ 
        return self(x)

    def get_field(self):
        r"""
            Return the finite field

        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3,sbox=[4,5,6,7,0,1,2,3])
            sage: S.get_field()
            Finite Field in a of size 2^3

            sage: S.get_field().fetch_int(7)
            a^2 + a + 1
        """
        return self._K

    def get_linear_functions(self):
        r"""
            Return [M1,M2,M3,V1,V2] from the last EA-equivalence ``F(x) = M1 * G( M2 * x + V2 ) + M3 * x + V1``

        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3)
            sage: S.generate_sbox(method="polynomial",G="x^3",T='EA')
            sage: M1,M2,M3,V1,V2 = S.get_linear_functions()
            sage: sbox = S.get_sbox()
            sage: S.generate_sbox(method="polynomial",G="x^3",T='EA',M1=M1,M2=M2,M3=M3,V1=V1,V2=V2)
            sage: S.get_sbox() == sbox
            True

            sage: S.generate_sbox(method="polynomial",G="x^3",T='EA',M1=M1,M2=M2,M3=M3,V1=V2,V2=V1)
            sage: S.get_sbox() == sbox
            False

            sage: S.generate_sbox(method="polynomial",G="x^3",T='EA',M1=M1,M2=M2)
            sage: S.get_sbox() == sbox
            False

            sage: [M1,M2] == S.get_linear_functions()[:2]
            True

            sage: S.generate_sbox(method="polynomial",G="x^3",T='EA')
            sage: S.get_sbox() == sbox
            False

        """
        return self._LFoEA

    def get_mg(self):
        r"""
            Return multiplicative generator of the finite field

        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3)
            sage: S.get_mg()^(2^3-1) == 1
            True
        """
        return self._K.multiplicative_generator()

    def get_modulus(self):
        r"""
            Return the irreducible polynomial of the field

        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3)
            sage: S.get_modulus()
            x^3 + x + 1
        """
        return self._K.modulus()

    def get_ring(self):
        r"""
            Return the ring

        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3)
            sage: S.get_ring()
            Univariate Polynomial Ring in x over Finite Field in a of size 2^3
        """
        return self._P

    def get_sbox(self):
        r"""
            Return SBox as integer array

        EXAMPLE::
            
            sage: S = Sbox(n=3,m=3,sbox=[0,1,5,7,2,3,4,6])
            sage: S.get_sbox() == [0,1,5,7,2,3,4,6]
            True
        """
        return self._S

    def get_polynomial(self):
        r"""
            Return vectorial boolean function correspoinding to the substitution

        EXAMPLE::

            sage: S = Sbox(n=3,m=3,sbox=[0, 1, 3, 4, 5, 6, 7, 2])
            sage: S.get_polynomial()
            'x^3'
        """
        if self._polynomial is None:
            self._polynomial = self.interpolation_polynomial()
        return self._polynomial

    def set_sbox(self,sbox=None):
        r"""
        Set S-box to predefined S-box
        
        INPUT::
        
            - ``sbox`` -- integer array defining transformation

        EXAMPLE::

            sage: S = Sbox(n=3,m=3)
            sage: S
            []

            sage: S.set_sbox([0, 1, 3, 4, 5, 6, 7, 2])
            sage: S
            [0, 1, 3, 4, 5, 6, 7, 2]
        """
        if sbox is None:
            raise TypeError("sbox is not defined")

        if len(sbox) != self._length:
            raise TypeError("The length ({0}) of given S-box differs from internal ({1})".format(len(sbox),self._length))

        if ZZ(max(sbox)).nbits() > self._m:
            raise TypeError("The bit length ({0}) of the maximum element of given S-box larger then internal 'm' ({1})".format(ZZ(max(sbox)).nbits(),self._m))

        self._polynomial = None
        self._S = sbox[:]
        self._LFoEA = [None,None,None,None,None]

    def l2m(self,L=None):
        r'''
        Convert a linear function ``L`` to a matrix ``M``
        
        INPUT::
        
            - ``L`` -- linear function to convert

        EXAMPLE::

            sage: S = Sbox(n=3,m=3)
            sage: S.l2m("g^5*x^8+g*x^4+x+g")
            [0 0 1]
            [1 1 1]
            [1 0 0]
        '''
        if L is None:
            raise TypeError("L is not defined")

        if isinstance(L,basestring):
            L = self.g2p(L)

        M=matrix(GF(2),self._n,self._m)

        for i in xrange(self._n):
            M.set_column(i,L.subs(self._K(ZZ(2^i).digits(2)))._vector_())

        return M

    def m2l(self,M=None,representation="generator"):
        r'''
        Convert a matrix ``M`` to a linear function ``L``

        INPUT::

            - ``M`` -- matrix to be converted
            - ``representation`` - a string, defines how the function ``L`` should be represented
                - 'generator' (default) -- the linear function will be given by it's generator
                - 'internal' -- the linear function will be given as internal representation

        EXAMPLE::

            sage: S = Sbox(n=3,m=3)
            sage: M = matrix(GF(2),3,[[0,0,1],[1,1,1],[1,0,0]])
            sage: S.m2l(M)
            'g^3*x^4 + g^4*x^2 + g^3*x'
        '''
        if M is None:
            raise TypeError("M is not defined")

        T = matrix(self._K,self._n)
        C = vector(self._K,self._n,[self._K(M.column(g)) for g in xrange(self._n)])

        for i in xrange(self._n):
            T.set_row(i,[self._K.fetch_int(2^i)^(2^g) for g in xrange(self._n)])

        T = (T.inverse()*C).list()

        L = sum([self._P("({0})*x^(2^{1})".format(T[g],g)) for g in xrange(len(T))])

        if representation == "internal":
            return L
        else:
            return self.p2g(pol=L)

    def p2g(self,pol=None,**kwargs):
        r"""
        Convert internal representation to polynomial with the primitive element.

        INPUT::
        
        - ``pol`` -- input polynomial
        - ``selftest`` -- if "True" (default) then function computes selftesting    

        EXAMPLE::

            sage: S=Sbox(n=3,m=3)
            sage: S.generate_sbox(method='random_substitution')
            sage: S.interpolation_polynomial() == S.p2g(S.g2p(S.interpolation_polynomial()))
            True

        """
        
        pol  = self._P(pol).mod(self._P("x^{0}+x".format(self._length)))
        selftest = kwargs.get('selftest',True)
        
        if pol is None:
            return "None"

        g_pol = ""
        for i in xrange(self._length-1,1,-1):
            t = pol[i].log(self._K.multiplicative_generator())
            if pol[i] != 0:
                if t == 1:
                    g_pol += "g*x^{0} + ".format(i)
                elif t == 0:
                    g_pol += "x^{0} + ".format(i)
                else:
                    g_pol += "g^{0}*x^{1} + ".format(t,i)

        t = pol[1].log(self._K.multiplicative_generator())
        if pol[1] != 0:
            if t == 1:
                g_pol += "g*x + "
            elif t == 0:
                g_pol += "x + "
            else:
                g_pol += "g^{0}*x + ".format(t)

        t = pol[0].log(self._K.multiplicative_generator())
        if pol[0] != 0:
            if t == 1:
                g_pol += "g"
            elif t == 0:
                g_pol += "1"
            else:
                g_pol += "g^{0}".format(t)
        else:
            if g_pol[-3:] == " + ":
                g_pol = g_pol[:-3]

        if g_pol == "":
            g_pol = "0"
            
        if selftest == True:
            if self._P(g_pol.replace("g","({0})".format(self._K.multiplicative_generator()))) == pol:
                return g_pol
            else:
                raise TypeError("You have found a bug!")
        else:
            return g_pol

    def random_permutation(self):
        r"""
        Generate random permutation. See  :meth:`SBox._random_permutation`.

        EXAMPLE::

            sage: S=Sbox(n=3,m=3)
            sage: S.random_permutation()
            sage: S # random
            [5, 3, 4, 6, 0, 1, 7, 2]

            sage: S.is_bijection()
            True
        """
        self.generate_sbox(method='random_permutation')

    def random_substitution(self):
        r"""
        Generate random permutation. See  :meth:`SBox._random_permutation`.

        EXAMPLE::

            sage: S=Sbox(n=3,m=3)
            sage: S.random_substitution()
            sage: S # random
            [0, 4, 1, 3, 1, 3, 7, 5]

        """
        self.generate_sbox(method='random_substitution')

    def Tr(self,x=None,n=None,m=None):
        r"""
        Return trace from ``F_{2^n}`` to ``F_{2^m}``

        INPUT:

            - ``x`` -- input value in ``F_{2^n}``
            - ``n`` -- degree of the field extension "from"
            - ``m`` -- degree of the field extension "to"

        EXAMPLE::

            sage: S=Sbox(n=4,m=4)
            sage: K=S.get_field()
            sage: S.Tr(K.fetch_int(13))
            1

            sage: S.Tr(K.fetch_int(13),m=2)
            a^2 + a
        """

        if m is None:
            m=1

        if n is None:
            n=x.parent().degree()

        return sum([x^(2^(i*m)) for i in xrange(n/m)])

    def Tr_pol(self,x=None,n=None,m=None):
        r"""
        Return trace polynomial from ``F_{2^n}`` to ``F_{2^m}``

        INPUT::

            - ``x`` -- input polynomial
            - ``n`` -- degree of the field extension "from"
            - ``m`` -- degree of the field extension "to"

        EXAMPLE::

            sage: S=Sbox(n=4,m=4)
            sage: K=S.get_field()
            sage: P=S.get_ring()
            sage: F='g*x^5 + g^2*x^2 + g^4'
            sage: G=S.g2p(F)
            sage: G
            a*x^5 + a^2*x^2 + a + 1

            sage: F == S.p2g(G)
            True

            sage: tr = S.Tr_pol(G,m=2)
            sage: tr
            (a^2 + 1)*x^8 + x^5 + a^2*x^2 + 1

            sage: S.Tr(G.subs(K.fetch_int(3)),m=2)
            a^2 + a + 1

            sage: tr.subs(K.fetch_int(3))
            a^2 + a + 1
        """
        
        if m is None:
            m=1

        if n is None:
            n=x.parent().base().degree()

        pol=sum([x^(2^(i*m)) for i in xrange(n/m)])
        pol = pol.mod(self._P("x^{0}+x".format(self._length)))

        return pol

    # Functions from CSbox.sage

    absolute_indicator=cr_absolute_indicator
    algebraic_immunity_sbox=cr_algebraic_immunity_sbox
    check_polynomial=cr_check_polynomial
    check_system=cr_check_system
    CI=cr_CI
    create_system=cr_create_system
    cycles=cr_cycles
    difference_distribution_matrix=cr_difference_distribution_matrix
    interpolation_polynomial=cr_interpolation_polynomial
    is_APN=cr_is_APN
    is_balanced=cr_is_balanced
    is_bijection=cr_is_bijection
    is_CCZ_equivalent=cr_is_CCZ_equivalent
    MDT=cr_MDT
    maximal_difference_probability=cr_maximal_difference_probability
    maximal_linear_bias=cr_maximal_linear_bias
    MLT=cr_MLT
    minimum_degree=cr_minimum_degree
    NL=cr_NL
    PC=cr_PC
    resilient=cr_resilient
    SAC=cr_SAC
    SSI=cr_SSI

    IsEquivalentToPermutation=cr_IsEquivalentToPermutation # experemental
    is_EA_equivalent=cr_is_EA_equivalent # experemental

    # Functions from GSbox.sage

    _APN6=gen_APN6
    _CCZ=gen_CCZ
    _dickson=gen_dickson
    _dobbertin=gen_dobbertin
    _EA=gen_EA
    _gold=gen_gold
    _inverse=gen_inverse
    _kasami=gen_kasami
    _niho=gen_niho
    _OP4=gen_OP4
    _random_permutation=gen_random_permutation
    _random_substitution=gen_random_substitution
    _welch=gen_welch

    _claude_matrix=gen_e_claude_matrix # experemental
    _AI=gen_e_AI # experemental
    _lilia=gen_e_lilia # experemental
    _lilia_base=gen_e_lilia_base # experemental
