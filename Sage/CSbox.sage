r"""
    Cryptanalysis functions

AUTHORS:

- Oleksandr Kazymyrov (2011-06-21): initial version

- Oleksandr Kazymyrov (2013-08-04): rewrite lots of functions

- Oleksandr Kazymyrov (2013-08-15): rewrite output of the "cycles" function

"""

#*****************************************************************************
#       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***************************************************************************** 

def cr_absolute_indicator(self):
    r"""
    Return the maximum value of the autocorrelation spectrum

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='polynomial',G="g*x^3+1")
        sage: S.absolute_indicator()
        8
    """

    return c_absolute_indicator(self._S,self._length)

def cr_algebraic_immunity_sbox(self):
    r"""
    Return algebraic immunity, that is the maximum degree of a system of equations that describes the S-box, and number of equations

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: sage: S.generate_sbox(method='polynomial',G="g*x^3+1")
        sage: S.algebraic_immunity_sbox()
        [2, 14]
    """

    return c_algebraic_immunity_sbox(self._S,self._length,self._n,self._m)

def cr_is_balanced(self):
    r"""
    Check S-box on balancedness

    EXAMPLE::

        sage: S=Sbox(n=3,m=3,sbox=[5,0,1,1,2,3,7,6])
        sage: S.is_balanced()
        False

        sage: S.generate_sbox(method='random_permutation')
        sage: S.is_balanced()
        True
    """

    return c_is_balanced(self._S,self._length,self._m)

def cr_check_polynomial(self):
    r"""
    Check on equality polynomial and lookup table representations

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='random_substitution')
        sage: S.check_polynomial()
        True
    """

    if self._polynomial is None:
        self.interpolation_polynomial(representation="internal")

    if [(self._polynomial.subs(self._K(ZZ(i).digits(2)))).integer_representation() for i in xrange(self._length)] == self._S:
        return True
    else:
        return False

def cr_check_system(self, system=None, degree=2):
    r"""
    Checks the system of equations describing the S-box

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='random_substitution')
        sage: S.check_system()
        True
    """
    P = BooleanPolynomialRing(self._n+self._m, ["x%d"%i for i in range(self._n)] + ["y%d"%i for i in range(self._m)])
    gens = [P("x%d"%(self._n-i-1)) for i in range(self._n)] + [P("y%d"%(self._m-i-1)) for i in range(self._m)]

    if system is None:
        system = self.create_system(degree=degree)

    if len(system) == 0:
        return "System doesn't exist"

    solutions = []
    for i in range(self._length):
        solutions.append( dict(zip(gens, list(reversed(ZZ(i).digits(base=2,padto=self._n))) + list(reversed(ZZ(self._S[i]).digits(base=2,padto=self._m))) )) )

    if any(f.subs(s) for f in system for s in solutions) == False:
        return True

    return False

def cr_CI(self):
    r"""
    Return the correlation immunity of the vectorial boolean function 

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='polynomial',G="x^3")
        sage: S.CI()
        0
    """
    return c_CI(self._S,self._length)

def cr_create_system(self, degree=2, groebner=False):
    r"""
    Create the system of equations describing the S-box

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='polynomial',G="x^3")
        sage: S.create_system()
        [x0*y1 + x2 + y0*y2,
         x0*y1 + x0*y2 + x1 + y0*y1,
         x0*y0 + x0*y2 + x0 + x1*y1 + y0*y1,
         x0*y2 + y0*y2 + y2,
         x0*y1 + y0*y1 + y1,
         x0*y0 + x0*y1 + x1*y1 + y0*y2 + y0 + y1*y2,
         x0*y1 + x1*x2 + y1*y2,
         x0*x2 + x0*y1,
         x0*y1 + x1*y1 + x2*y2 + y0*y1 + y0*y2,
         x1*y1 + x2*y1 + y0*y1 + y1*y2,
         x2*y0 + y0*y2,
         x0*x1 + x0*y1 + x0*y2,
         x0*y2 + x1*y2 + y1*y2,
         x1*y0 + y0*y1]
    """
    #t1=time()

    P = BooleanPolynomialRing(self._n+self._m, ["x%d"%i for i in range(self._n)] + ["y%d"%i for i in range(self._m)])
    X = [P("x%d"%(self._n-i-1)) for i in range(self._n)]
    Y = [P("y%d"%(self._m-i-1)) for i in range(self._m)]

    gens = X+Y

    bits = []
    for i in xrange(self._length):
        bits.append( list(reversed(ZZ(i).digits(base=2,padto=self._n))) + list(reversed(ZZ(self._S[i]).digits(base=2,padto=self._m))) )

    nrows = self._length
    ncols = sum(binomial(self._m+self._n,i) for i in range(0,degree+1))

    A = Matrix(GF(2), nrows , ncols)

    exponents = []
    for d in xrange(degree+1):
        exponents += IntegerVectors(d, max_length=self._n+self._m, min_length=self._n+self._m, min_part=0, max_part=1).list()

    col = 0
    variables = []
    for exponent in exponents:
        variables.append( mul([gens[i]**exponent[i] for i in range(len(exponent))]))
        for row in xrange(self._length):
            A[row,col] = mul([bits[row][i] for i in range(len(exponent)) if exponent[i]])
        col +=1

    system=A.right_kernel()
    system=system.matrix()

    gens=[]
    length=len(variables)
    for j in xrange(len(system.rows())):
        gens.append(sum(variables[i]*system[j][i] for i in xrange(length)) )

    return gens

    #if not groebner:
    #    return gens

    #FI = set(FieldIdeal(P).gens())
    #I = Ideal(gens + list(FI))
    #gb = I.groebner_basis()

    #gens = []
    #for f in gb:
    #    if f not in FI: # filter out field equations
    #        gens.append(f)

    #return gens

def cr_cycles(self,**kargs):
    r"""
    Return cycle structure of the substitution

    INPUT::
        
    - ``graph`` -- if ``True`` save cycles in ``dot`` format, which can be imported by gephi (see http://www.gephi.org)
    - ``file`` -- path to the file with a graph (default is ``cycles.dot``)
    - ``format`` -- string, how to represent the result, e.g.,
        - 'chains' (default) -- all chains in full representation
        - 'chains_compact' -- only chains in the form [ [start point_1,length_1], [start point_2, length_2], ... ]
        - 'cycles' -- only cycles
        - 'cycles_compact' -- only cycles in the form [ [start point_1,length_1], [start point_2, length_2], ... ]
        - 'partition' -- [number of chains, number of pure chains, number of cycles, number of pure cycles]
        - 'pchains' -- only pure chains, i.e. chains without cycles
        - 'pchains_compact' -- only pure chains in the form [ [start point_1,length_1], [start point_2, length_2], ... ]
        - 'pcycles' -- only pure cycles, i.e. the lenght of a cycle equals the lenght of a chain
        - 'pcycles_compact' -- only pure cycles in the form [ [start point_1,length_1], [start point_2, length_2], ... ]

    EXAMPLE::

        An example of a substitution with all outputs:

::

        sage: S=Sbox(n=3,m=6,sbox=[0,3,1,2,32,55,1,5])
        sage: S.cycles()
        [[0, 0], [4, 32], [7, 5, 55], [6, 1, 3, 2, 1]]

        sage: S.cycles(format='chains') == S.cycles()
        True

        sage: S.cycles(format='chains_compact')
        [[0, 1], [1, 3]]

        sage: S.cycles(format='cycles')
        [[0, 0], [1, 3, 2, 1]]

        sage: S.cycles(format='cycles_compact')
        [[0, 1], [1, 3]]

        sage: S.cycles(format='partition')
        [4, 2, 1]

        sage: S.cycles(format='pchains')
        [[4, 32], [7, 5, 55]]

        sage: S.cycles(format='pchains_compact')
        [[4, 2], [7, 3]]

        sage: S.cycles(format='pcycles')
        [[0, 0]]

        sage: S.cycles(format='pcycles_compact')
        [[0, 1]]

        An example of a substitution with one cycle:

::
        sage: S=Sbox(n=3,m=6,sbox=[3, 0, 6, 4, 1, 0, 0, 1])
        sage: S.cycles()
        [[5, 0, 3, 4, 1, 0], [7, 1, 0, 3, 4, 1], [2, 6, 0, 3, 4, 1, 0]]

        sage: S.cycles(format='chains_compact')
        [[5, 5], [7, 5], [2, 6]]

        sage: S.cycles(format='cycles')
        [[0, 3, 4, 1, 0]]

        sage: S.cycles(format='cycles_compact')
        [[0, 4]]

        sage: S.cycles(format='partition')
        [3, 0, 0]

        sage: S.cycles(format='pchains')
        []

        sage: S.cycles(format='pchains_compact')
        []

        sage: S.cycles(format='pcycles')
        []

        sage: S.cycles(format='pcycles_compact')
        []
    """
    graph=kargs.get('graph',False)
    format=kargs.get('format',"chains")

    cycles = cpp_cycles(self._S,self._length)

    if graph == True:
        e = []
        for c in cycles:
            for a,b in [(c[i], c[i+1]) for i in xrange(len(c)-1)]:
                e.append([a,b])

        e = [list(i) for i in set([tuple(i) for i in e])]

        g = DiGraph({}, loops=True, multiedges=True)
        for i in xrange(len(e)):
            g.add_edge(e[i][0], e[i][1])

        g.graphviz_to_file_named(kargs.get('file','cycles.dot'))

    output = []

    if format == "chains":

        return cycles

    elif format == "chains_compact":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k != len(c)-1:
                output.append([c[0],len(c)-1])
            else:
                output.append([c[0],len(c)])

    elif format == "cycles":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k != len(c)-1:
                # new cycle will start with the lowest number
                tmp = c[k:]

                m = tmp.index(min(tmp)) # find the minimum numer
                tmp = c[k:] + c[k+1:] # double cycle
                tmp = tmp[m:m+len(c[k:])] # copy the new cycle

                if not tmp in output:
                    output.append(tmp)

    elif format == "cycles_compact":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k != len(c)-1:
                # new cycle will start with the lowest number
                tmp = c[k:]

                m = tmp.index(min(tmp)) # find the minimum numer
                tmp = c[k:] + c[k+1:] # double cycle
                tmp = tmp[m:m+len(c[k:])] # copy the new cycle

                tmp = [tmp[0],len(tmp)-1]

                if not tmp in output:
                    output.append(tmp)

    elif format == "partition":

        output = [0,0,0,0]
        unique_cycles = []

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k != len(c)-1:
                # new cycle will start with the lowest number
                tmp = c[k:]

                m = tmp.index(min(tmp)) # find the minimum numer
                tmp = c[k:] + c[k+1:] # double cycle
                tmp = tmp[m:m+len(c[k:])] # copy the new cycle

                if not tmp in unique_cycles:

                    if c[0] == c[-1]:
                        output[3] += 1
                    else:
                        output[2] += 1
                            
                    unique_cycles.append(tmp)
            else:
                output[1] += 1

            output[0] += 1

    elif format == "pchains":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k == len(c)-1:
                output.append(c)

    elif format == "pchains_compact":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k == len(c)-1:
                output.append([c[0],len(c)])

    elif format == "pcycles":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k != len(c)-1:
                if c[0] == c[-1]:
                    # new cycle will start with the lowest number
                    tmp = c[k:]

                    m = tmp.index(min(tmp)) # find the minimum numer
                    tmp = c[k:] + c[k+1:] # double cycle
                    tmp = tmp[m:m+len(c[k:])] # copy the new cycle

                    if not tmp in output:
                        output.append(tmp)

    elif format == "pcycles_compact":

        for i,c in enumerate(cycles):
            k = c.index(c[-1])

            if k != len(c)-1:
                if c[0] == c[-1]:
                    # new cycle will start with the lowest number
                    tmp = c[k:]

                    m = tmp.index(min(tmp)) # find the minimum numer
                    tmp = c[k:] + c[k+1:] # double cycle
                    tmp = tmp[m:m+len(c[k:])] # copy the new cycle

                    tmp = [tmp[0],len(tmp)-1]

                    if not tmp in output:
                        output.append(tmp)

    else:
        raise TypeError("Unsupported format '{0}'".format(format))

    return output

def cr_difference_distribution_matrix(self):
    r"""
    Return the minimum degree of substitution.

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='polynomial',G="x^3")
        sage: S.difference_distribution_matrix()
        [8 0 0 0 0 0 0 0]
        [0 2 0 2 0 2 0 2]
        [0 0 2 2 2 2 0 0]
        [0 2 2 0 2 0 0 2]
        [0 0 0 0 2 2 2 2]
        [0 2 0 2 2 0 2 0]
        [0 0 2 2 0 0 2 2]
        [0 2 2 0 0 2 2 0]

        sage: S.random_substitution()
        sage: S.difference_distribution_matrix() # random
        [8 0 0 0 0 0 0 0]
        [0 4 2 0 0 0 2 0]
        [4 2 0 0 0 2 0 0]
        [0 4 0 2 0 0 0 2]
        [2 0 0 0 0 2 2 2]
        [0 2 0 0 2 0 2 2]
        [2 0 0 0 0 2 2 2]
        [0 2 0 0 2 0 2 2]

        sage: S=Sbox(n=4,m=3)
        sage: S.random_substitution()
        sage: S.difference_distribution_matrix() # random
        [16  0  0  0  0  0  0  0]
        [ 2  0  2  4  0  2  6  0]
        [ 4  0  2  2  0  4  4  0]
        [ 4  0  0  0  2  2  8  0]
        [ 4  0  2  2  0  4  4  0]
        [ 6  0  0  2  2  0  6  0]
        [ 8  0  0  0  2  2  4  0]
        [ 0  0  0  4  2  2  8  0]
        [ 4  2  0  6  0  2  2  0]
        [ 4  0  2  2  0  8  0  0]
        [ 2  0  0  6  2  4  2  0]
        [ 2  0  0  2  0  6  4  2]
        [ 4  2  0  6  0  2  2  0]
        [ 4  0  0  4  0  4  2  2]
        [ 0  2  0  6  0  2  6  0]
        [ 2  0  0  2  2  8  2  0]

        sage: S=Sbox(n=3,m=4)
        sage: S.random_substitution()
        sage: S.difference_distribution_matrix() # random
        [8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 2 2 0 0 0 0 0 0 0 0 2 0 0 2]
        [2 0 0 2 0 0 0 0 0 0 0 0 2 2 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 2 0 2 4]
        [2 0 2 0 0 0 0 0 0 0 0 0 0 0 0 4]
        [2 0 0 2 0 0 0 0 0 0 0 0 2 2 0 0]
        [0 0 2 2 0 0 0 0 0 0 0 0 2 0 0 2]
        [4 0 0 0 0 0 0 0 0 0 0 0 2 0 2 0]

    NOTE:: 

        This functions is very slow and can be used only for small substitutions
    """
    nrows = 1<<self._n
    ncols = 1<<self._m

    A = Matrix(ZZ, nrows, ncols)

    for x1 in xrange(nrows):
        for x2 in xrange(nrows):
            A[ x1^^x2 , self._S[x1]^^self._S[x2]] += 1

    return A

def cr_interpolation_polynomial(self, representation="generator"):
    r"""
    Return the interpolation polynomial of the substitution.

    INPUT:

        - ``representation`` -- string, how to represent result, e.g.,
            - ``generator`` (default) -- return polynomial with generators
            - ``internal`` -- return polynomial in iternal representation

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='polynomial',G='g^6*x^5+x^2+x+g^7')
        sage: S.interpolation_polynomial()
        'g^6*x^5 + x^2 + x + 1'

    """
    if self._polynomial is not None:
        if representation == "internal":
            return self._polynomial
        else:
            return self.p2g(pol=self._polynomial)

    l = []
    for i in xrange(self._length):
        l.append( (self._K(ZZ(i).digits(2)), self._K(ZZ(self._S[i]).digits(2))) )

    self._polynomial = self._P.lagrange_polynomial(l)

    if representation == "internal":
        return self._polynomial
    else:
        return self.p2g(pol=self._polynomial)

def cr_is_bijection(self):
    r"""
    Check substitution on bijection

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='random_substitution')
        sage: S.is_bijection()
        False

        sage: S.generate_sbox(method='random_permutation')
        sage: S.is_bijection()
        True

        sage: S=Sbox(n=6,m=3,sbox=[0, 1, 3, 4, 5, 6, 7, 2]*8)
        sage: S.is_bijection()
        False

        sage: S.is_balanced()
        True
    """
    if self._n == self._m:
        return len(set(self._S)) == (self._length)
    else:
        return False

def cr_is_APN(self, mode="c"):
    r"""
    Return ``True`` if the substitution is APN (2-uniform)

    INPUT:
    
    - ``mode`` -- string, which algorithm is used for calculation, e.g.,
        - ``c`` (default) -- "C" program 
        - ``sage`` -- Sage implementation 

    EXAMPLE::

        sage: S=Sbox(n=6,m=6)
        sage: S.generate_sbox(method='APN6')
        sage: S.is_APN()
        True

        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.is_APN()
        False
    """

    if mode == "sage":
        T = zero_vector(ZZ,self._length)

        for d in xrange(1,self._length):
            for x in xrange(self._length):
                T[self._S[x]^^self._S[x^^d]] += 1
            if len([i for i in T if i > 2]) != 0:
                return False
            T=zero_vector(ZZ,self._length)
        return True
    else:
        if c_is_APN(self._S,self._length,max(1<<self._m,self._length)) == True:
            return True
        return False

def cr_is_CCZ_equivalent(self,F=None,G=None):
    r"""
    Return ``True`` if ``F`` and ``G`` are CCZ-equivalent. Testing is based on
    the isomorphism of codes.

    - ``F,G`` -- two functions for comparison

    EXAMPLE::

        sage: S=Sbox(n=3,m=3)
        sage: S.is_CCZ_equivalent('x^3','x^4')
        False

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='random_substitution')
        sage: F=S.interpolation_polynomial()
        sage: S.generate_sbox(method='polynomial',G=F,T='CCZ')
        sage: G=S.interpolation_polynomial()
        sage: S.is_CCZ_equivalent(F,G)
        True

    """

    if F is None:
        raise TypeError("Function 'F' must be presented")

    if G is None:
        raise TypeError("Function 'G' must be presented")

    polF=self.g2p(F)
    polG=self.g2p(G)

    M=Matrix(GF(2), self._length, self._n + self._m + 1, [[1]+ZZ(i).digits(2,padto=self._n)+vector(polF.subs(self._K.fetch_int(i))).list() for i in xrange(self._length)])
    M=M.transpose()
    CF=LinearCode(M)
    
    M=Matrix(GF(2), self._length, self._n + self._m + 1, [[1]+ZZ(i).digits(2,padto=self._n)+vector(polG.subs(self._K.fetch_int(i))).list() for i in xrange(self._length)])
    M=M.transpose()
    CG=LinearCode(M)

    return CF.is_permutation_equivalent(CG)

def cr_fixed_points(self):
    r"""
    Return the list of fixed points

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.random_substitution()
        sage: S.fixed_points() # random
        []

        sage: S=Sbox(n=3,m=3)
        sage: S.set_sbox([0,5,1,3,4,2,6,7])
        sage: S.fixed_points()
        [0,3,4,6,7]
    """

    return cpp_fixed_points(self._S,self._length)

def cr_MDT(self,difference='XORXOR'):
    r"""
    Return the maximum value of a differential table

    INPUT::
    
        - ``difference`` –- string, describes a key addition method:
                            ``XORXOR`` (default) --  compute difference of XOR->XOR 
                            ``XORADD``           --  compute difference of XOR->ADD 
                            ``ADDXOR``           --  compute difference of ADD->XOR 
                            ``ADDADD``           --  compute difference of ADD->ADD

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.MDT()
        2

        sage: S.MDT("ADDXOR")
        4

        sage: S.MDT("XORADD")
        4

        sage: S.MDT("ADDADD")
        4

        sage: S=Sbox(n=8,m=8)
        sage: S.generate_sbox(method='inverse')
        sage: S.MDT()
        4

        sage: S.MDT("ADDXOR")
        6

        sage: S.MDT("XORADD")
        6

        sage: S.MDT("ADDADD")
        6

        sage: M1=matrix(GF(2),8,[[0,0,0,0,1,1,0,1],[1,0,0,0,1,1,1,1],[0,1,1,0,0,1,0,1],[1,0,1,0,1,1,0,1],[1,0,0,1,0,1,0,0],[1,1,1,0,1,0,1,1],[0,1,1,0,1,1,0,1],[0,0,0,0,1,1,1,1]])
        sage: M2=matrix(GF(2),8,[[0,1,1,0,0,0,0,1],[1,0,0,1,0,1,0,1],[0,0,0,0,1,1,0,0],[0,0,1,1,0,0,0,0],[0,1,1,0,1,0,1,0],[0,0,0,1,1,0,0,0],[0,1,0,0,1,1,0,1],[1,1,1,0,0,1,1,0]])
        sage: V1=vector(GF(2),8,[0, 1, 0, 0, 0, 1, 1, 1])
        sage: V2=vector(GF(2),8,[0, 0, 1, 1, 0, 0, 0, 0])
        sage: S.generate_sbox(method='inverse',T="A",M1=M1,M2=M2,V1=V1,V2=V2)
        sage: S.MDT()
        4

        sage: S.MDT("ADDXOR")
        7

        sage: S.MDT("XORADD")
        6

        sage: S.MDT("ADDADD")
        7
    """

    if difference not in ['XORXOR','XORADD','ADDXOR','ADDADD']:
        raise NotImplementedError("'{0}'".format(difference))

    return c_MDT(self._S,self._n,self._m,difference)

def cr_maximal_difference_probability(self,difference='XORXOR'):
    r"""
    Return the maximum value difference probability

    INPUT::
    
        - ``difference`` –- string (``XORXOR`` (default), ``XORADD``, ``ADDXOR`` and ``ADDADD``) describes a key addition method

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.maximal_difference_probability()
        0.0312500000000000

        sage: S.maximal_difference_probability('XORADD')
        0.0625000000000000

        sage: S.maximal_difference_probability('ADDXOR')
        0.0625000000000000

        sage: S.maximal_difference_probability('ADDADD')
        0.0625000000000000
    """

    return self.MDT(difference)/(2.0<<self._n)

def cr_maximal_linear_bias(self,addition='XOR'):
    r"""
    Return the maximum linear bias

    INPUT::
    
    - addition –- string, describes a key addition method, e.g.,
        - ``'XOR'`` --  compute difference for XOR->XOR (default)

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.maximal_linear_bias()
        0.0625000000000000

    TODO::

        Implement ``XORADD``, ``ADDXOR`` and ``ADD``.
    """

    if addition != 'XOR':
        raise NotImplementedError("'{0}'".format(addition))

    return self.MLT(addition)/(2.0<<self._n)

def cr_minimum_degree(self):
    r"""
    Return the minimum degree of the substitution

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.minimum_degree()
        2
    """
    return c_minimum_degree(self._S,self._n,self._m)

def cr_MLT(self,addition='XOR'):
    r"""
    Return the maximum value of linear approximation table.

    INPUT::
    
    - addition –- string, describes a key addition method, e.g.,
        - ``'XOR'`` --  compute difference for XOR->XOR (default)

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.MLT()
        4

    TODO::

        Implement ``XORADD``, ``ADDXOR`` and ``ADD``.
    """

    if addition != 'XOR':
        raise NotImplementedError("'{0}'".format(addition))

    return c_MLT(self._S,self._n,self._m)

def cr_NL(self):
    r"""
    Return the nonlinearity of the substitution

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.NL()
        12
    """
    return (self._length>>1)-c_MLT(self._S,self._n,self._m)

def cr_PC(self):
    r"""
    Return the value of propogation criterion

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.PC()
        0
    """
    return c_PC(self._S,self._length)

def cr_resilient(self):
    r"""
    Return resiliency of the substitution.

    If function is not balanced then return ``False``.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.resilient()
        0 

        sage: S=Sbox(n=3,m=3,sbox=[0,1,2,3,4,5,6,0])
        sage: S.resilient()
        False

    """

    if c_is_balanced(self._S,self._length,self._m):
        return c_CI(self._S,self._length)

    return False

def cr_SAC(self):
    r"""
    Return ``True`` if the substitution satisfies strict avalanche criterion.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.SAC()
        False
    """

    if c_PC(self._S,self._length):
        return True
    else:
        return False

def cr_SSI(self):
    r"""
    Return the sum-of-squares indicator.

    EXAMPLE::

        sage: S=Sbox(n=5,m=5)
        sage: S.generate_sbox(method='polynomial',G='x^10')
        sage: S.SSI()
        2048
    """
    return c_SSI(self._S,self._length)
