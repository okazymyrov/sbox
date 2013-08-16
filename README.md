# Overview

A substitution box or S-box is one of the basic components of symmetric key cryptography. In general, an S-box takes _n_ input bits and transforms them into _m_ output bits. This is called an _n x m_ S-box and is often implemented as a lookup table. 

Another representation is by vectorial Boolean function, which is defined by a Polynomial Ring in _x_ over Finite Field in _a_ of size 2<sup>n</sup> with primitive element _g_.

Note that by default bits are interpreted in little endian notation, e.g.

> 6<sub>10</sub> = 110<sub>2</sub>

> 10<sub>10</sub> = 1010<sub>2</sub> = A<sub>16</sub>

# Authors

- Oleksandr Kazymyrov
- Maksim Storetvedt
- Anna Maria Eilertsen

# How to use

1. In Sage. Go to the ''Sage'' directory and type ``load ./Sbox.sage``.
2. Use externaly (see examples in ''Main.sage'' ).

# Examples

* Create an example of a _3 x 3_ substitution

		sage: S = Sbox(n=3,m=3,sbox=[1,2,3,0,7,6,5,4])
		sage: S.get_sbox()
		[1, 2, 3, 0, 7, 6, 5, 4]
		
		sage: F=S.interpolation_polynomial()
		sage: F
		'g^2*x^6 + g^3*x^5 + g^4*x^4 + g^5*x^3 + g^6*x^2 + 1'
		
		sage: S.generate_sbox(method='polynomial',G=F)
		sage: S.get_sbox() == [1, 2, 3, 0, 7, 6, 5, 4]
		True

* Generate APN permutation for _n = 6_ [1]

        sage: S=Sbox(n=6,m=6)
        sage: S.generate_sbox(method='APN6')
        sage: S.is_APN()
        True
        sage: S.MDT()
        2
        sage: S.is_bijection()
        True

* Apply CCZ and EA equivalence

        sage: S=Sbox(n=3,m=3)
        sage: S.generate_sbox(method='polynomial',G='x^3')
        sage: S
        [0, 1, 3, 4, 5, 6, 7, 2]
        
        sage: S.generate_sbox(method='polynomial',G='x^3',T='CCZ')
        sage: S # random
        [0, 6, 2, 3, 6, 3, 6, 4]
        sage: S.is_APN()
        True
        
        sage: S.generate_sbox(method='polynomial',G='x^3',T='EA')
        sage: S # random
        [7, 0, 7, 5, 0, 6, 6, 5]

        sage: S.is_APN()
        True
        		
# References

1. [J. F. Dillon, APN Polynomials: An Update. Fq9, International Conference on Finite Fields
and their Applications, University College Dublin. July 2009.](http://mathsci.ucd.ie/~gmg/Fq9Talks/Dillon.pdf)
2. [Sage](http://sagemath.org/)