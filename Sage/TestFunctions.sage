r"""
    Examples of using the class "Sbox"

AUTHORS:

- Oleksandr Kazymyrov (2011-06-21): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013 Oleksandr Kazymyrov <okazymyrov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***************************************************************************** 

load("./Data.sage")
load("./Sbox.sage")

def test_temp(**kargs):
    bits=kargs.get('bits',3)

    S=Sbox(n=bits,m=bits)

    S.generate_sbox(method='APN6')

    print "Sbox\t\t\t\t: {0}".format(S.get_sbox())
    print ""
    print "Characteristics of component functions"
    print "Balanced\t\t\t: {0}".format(S.is_balanced())
    print "Nonlinearity\t\t\t: {0}".format(S.NL())
    print "Absolute indicator\t\t: {0}".format(S.absolute_indicator())
    print "Propagation criterion\t\t: {0}".format(S.PC())
    print "Correlation immunity\t\t: {0}".format(S.CI())
    print "Sum-of-squares indicator\t: {0}".format(S.SSI())
    print "Minimum degree\t\t\t: {0}".format(S.minimum_degree())
    print "Resilient\t\t\t: {0}".format(S.resilient())
    print "SAC\t\t\t\t: {0}".format(S.SAC())
    print ""
    print "Characteristics of the substitution"
    print "Bijection\t\t\t: {0}".format(S.is_bijection())
    print "Interpolation polynomial\t: {0}".format(S.interpolation_polynomial())
    print "Check interpolation polynomial\t: {0}".format(S.check_polynomial())
    print "Multiplicative generator\t: {0}".format(S.get_mg())
    print "Modulus\t\t\t\t: {0}".format(S.get_modulus())
    print "Maximum of diff table\t\t: {0}".format(S.MDT())
    print "Maximum of lin table\t\t: {0}".format(S.MLT())
    print "Cycles\t\t\t\t: {0}".format(S.cycles())
    ret=S.algebraic_immunity_sbox()
    print "Algebraic immunity\t\t: degree={0} equations={1}".format(ret[0],ret[1])
    print "Check system\t\t\t: {0}".format(S.check_system(degree=ret[0]))
    print "~"*40
    print ""

def test_all_functions(**kargs):
    bits=kargs.get('bits',3)

    Sboxes = [
        AES_sbox
    ]

    S=Sbox(n=bits,m=bits)

    for sbox in Sboxes:
        S.set_sbox(sbox)
        print "Sbox\t\t\t\t:= {0}".format(S.get_sbox())
        print ""
        print "Characteristics of component functions"
        print "Balanced\t\t\t: {0}".format(S.is_balanced())
        print "Nonlinearity\t\t\t: {0}".format(S.NL())
        print "Absolute indicator\t\t: {0}".format(S.absolute_indicator())
        print "Propagation criterion\t\t: {0}".format(S.PC())
        print "Correlation immunity\t\t: {0}".format(S.CI())
        print "Sum-of-squares indicator\t: {0}".format(S.SSI())
        print "Minimum degree\t\t\t: {0}".format(S.minimum_degree())
        print "Resilient\t\t\t: {0}".format(S.resilient())
        print "SAC\t\t\t\t: {0}".format(S.SAC())
        print ""
        print "Characteristics of the substitution"
        print "Bijection\t\t\t: {0}".format(S.is_bijection())
        print "Interpolation polynomial\t: {0}".format(S.interpolation_polynomial())
        print "Check interpolation polynomial\t: {0}".format(S.check_polynomial())
        print "Multiplicative generator\t: {0}".format(S.get_mg())
        print "Modulus\t\t\t\t: {0}".format(S.get_modulus())
        print "Maximum of diff table\t\t: {0}".format(S.MDT())
        print "Maximum of lin table\t\t: {0}".format(S.MLT())
        print "Cycles\t\t\t\t: {0}".format(S.cycles())
        ret=S.algebraic_immunity_sbox()
        print "Algebraic immunity\t\t: degree={0} equations={1}".format(ret[0],ret[1])
        print "Check system\t\t\t: {0}".format(S.check_system(degree=ret[0]))
        print "~"*40
        print ""
                
    return 0

