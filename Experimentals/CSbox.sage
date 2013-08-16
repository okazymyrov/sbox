def cr_algebraic_immunity_bf(self):
    ai=0
    ret=self._n

    B = BooleanFunction(self._n)
    for j in xrange(1,1<<self._m):
        for i in xrange(1<<self._m
        ):
            B[i] = (self(i)&j).popcount()
        ai = B.algebraic_immunity()
        print "j={0}: ai={1}".format(j,ai)
        if ai < ret:
            ret=ai

	print "done"
	return ret

def cr_is_CCZ_equivalent(self,F=None,G=None,functions=False):
	r'''
		>>> Function has beta status <<<
	'''
	def check_pCCZ(sboxF,sboxG):
		r'''
			Return functions if esboxF_b(x) = esboxG_b(M * x), where esboxF_b is boolean function from F_2^{2n} to F_2.
		'''

		if c_is_CCZ_equivalent(sboxF,sboxG,self._length) != 0:
			return None
		else:
			print "do something..."
	
	if F is None:
		raise TypeError("Function 'F' must be presented")

	if G is None:
		raise TypeError("Function 'G' must be presented")

	sboxF = range(self._length)
	sboxG = range(self._length)
	
	sboxF = [F.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]
	sboxG = [G.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]

	M1 = None
	V1 = None
	M2 = None
	V2 = None
	M3 = None

	#for v2 in xrange(self._length):
	for v2 in xrange(1):
		sbox=sboxF[:]
		sbox=[sbox[g^^v2] for g in xrange(self._length)]
		ret=check_pCCZ(sbox,sboxG)

		if ret != None:
			M1=ret[0]
			V1=ret[1]
			V2=vector(GF(2),ZZ(v2).digits(base=2,padto=self._n))
			if functions == True:
				return [M1,V1,M2,V2,M3]
			else:
				return True

	return False

def cr_IsEquivalentToPermutation(self,**kwargs):
	r'''
		Is function 'F' equivalent to APN permutation?

		The algorithm is discribed in
		1. http://igitur-archive.library.uu.nl/student-theses/2012-0928-200735/MasterThesis-Eric_Cornet_29-aug-2012.pdf
		2. http://goo.gl/5pMHU
		
		Remarks:
			Is it necessary to use 'save' function? 
	'''
	F = kwargs.get('F',None)
	debug = kwargs.get('debug',False)
	save  = kwargs.get('save',False)
	retL  = kwargs.get('L',[])
	L 	  = kwargs.get('L1',matrix(GF(2),self._n,0))

	if F is None:
		raise TypeError("Function 'F' must be presented")

	if not isinstance(retL,list):
		raise TypeError("Ls must be in list")

	Sigma = [[] for _ in xrange(self._m+self._n)]

	self.generate_sbox(method='polynomial',G=F)

	if self._S[0] != 0:
		self._S = [self._S[0]^^self._S[g] for g in xrange(self._length)]

	for x in xrange(self._length):
		if debug:
			sys.stdout.write("\r[%-100s] %d%%" % ('='*(int(x*100/self._length)), (x*100/self._length).n(2)))
			sys.stdout.flush()
		for y in xrange(x+1,self._length):
			c = matrix(GF(2),self._n+self._m,1,flatten([ZZ(x^^y).digits(2,padto=self._n),ZZ(self._S[x]^^self._S[y]).digits(2,padto=self._m)]))
			ic = ZZ(c.list(),2).nbits()
			Sigma[ic-1].append(c[:ic,:])

	if debug:
		sys.stdout.write("\r[%-100s] %d%%\n" % ('='*100, 100))

	M=matrix(GF(2),self._m+self._n,self._m+self._n)
	z = zero_matrix(GF(2),self._n,1)	
	stop = 0
	
	if L != matrix(GF(2),self._n,0):
		j = L.ncols()
		indexes = [ZZ(g.list(),2)+1 for g in L.columns()] + [0 for _ in xrange(self._n+self._m-j)]
		retL = [g.echelon_form() for g in retL]
	else:
		j = 0
		indexes = [0 for _ in xrange(self._n+self._m)]

	while True:
		l = indexes[j]

		if L.ncols() == floor(self._n*1.4):
			print "check point = {0} ({1})".format([ZZ(g.list(),2) for g in L.columns()],len(retL))
		
		if L.ncols() <= j:
			L = matrix(GF(2),self._n,j+1,flatten([g.list()+[0] for g in L.rows()]))
		else:
			L = L[:,:j+1]

		while True:
			if j < self._n:
				if l == ((1<<j)+1):
					indexes[j] = 0
					j -= 1
					break
			else:
				if l == self._length:
					indexes[j] = 0
					j -= 1
					break

			L.set_column(j,vector(GF(2),ZZ(l).digits(2,padto=self._n)))

			if L.echelon_form() != L:
				l += 1
				continue

			next = 0

			if (j == (self._n+self._m-1)) and (L.rank() != self._n):
					next = 1
			else:
				for c in Sigma[j]:
					if L*c == z:
						next = 1
						break

			if next == 1:
				l += 1
			else:
				if j == self._n+self._m-1:
					eL = L.echelon_form()
					if not eL in retL:
						if debug:
							print "L = {0} ({1})".format([ZZ(g.list(),2) for g in eL.columns()],len(retL))
						for rL in retL:
							M.set_block(0,0,rL)
							M.set_block(self._m,0,eL)

							if not M.is_singular():
								stop = 1
								j = -1
								break
						if j == -1:
							break
						retL.append(copy(eL))
					l += 1
				else:
					indexes[j] = l + 1
					j += 1
					break
		if j == -1:
			break

	if stop == 0:
		return None
	else:
		return M

def cr_is_EA_equivalent(self,F=None,G=None,functions=False):
	r'''
		>>> Function has beta status! <<<
	'''
	def check_pEA1(sboxF,sboxG):
		r'''
			Return functions M and V if sboxF(x) = M * sboxG(x) + V
		'''
		M = matrix(GF(2),nrows=self._m,ncols=self._n)

		sboxFt=sboxF[:]
		sboxGt=sboxG[:]

		if sboxFt[0] != 0:
			sboxFt = [ g^^sboxFt[0] for g in sboxFt ]

		if sboxGt[0] != 0:
			sboxGt = [ g^^sboxGt[0] for g in sboxGt ]

		sboxGt_1=[ i/i-2 for i in xrange(1,self._length+1)]

		for i in xrange(len(sboxGt)):
			sboxGt_1[sboxGt[i]]=i

		V = vector(GF(2),ZZ(sboxF[0]).digits(base=2,padto=self._n))

		for i in xrange(self._n):
			x=sboxGt_1[1<<i]
			M.set_column(i,ZZ(sboxFt[x]).digits(base=2,padto=self._n))

		Ma = M*vector(GF(2),ZZ(sboxG[0]).digits(base=2,padto=self._n))		
		V = vector(GF(2),[ ZZ(Ma[g]) ^^ ZZ(V[g]) for g in xrange(self._n) ])

		sbox = range(self._length)

		V=ZZ(V.list(),2)

		for i in xrange(self._length):
			sbox[i]=sboxG[i]

			sbox[i]=vector(GF(2),ZZ(sbox[i]).digits(base=2,padto=self._n))
			sbox[i]=M*sbox[i]
			sbox[i]=ZZ(sbox[i].list(),2)
			
			sbox[i]= sbox[i]^^V

		V=vector(V.digits(2,padto=self._n))

		#print "~"*10
		#print "M:\n{0}".format(M)
		#print "V:\n{0}".format(V)
		#print "~"*10

		#print "sboxF  = {0}".format(sboxF)
		#print "sbox   = {0}".format(sbox)

		if sboxF == sbox:
			return [M,V]
		else:
			return None

	def check_pEA2(sboxF,sboxG):
		r'''
			Return functions M and V if sboxF(x) = sboxG(M * x + V)
		'''
		M = matrix(GF(2),nrows=self._m,ncols=self._n)

		sboxFt=sboxF[:]
		sboxGt=sboxG[:]

		sboxGt_1=[ -1 for i in xrange(1,self._length+1)]

		for i in xrange(len(sboxGt)):
			sboxGt_1[sboxGt[i]]=i

		if -1 in sboxGt_1:
			return None

		for i in xrange(self._length):
			sboxFt[i]=sboxGt_1[sboxF[i]]

		V = vector(GF(2),ZZ(sboxFt[0]).digits(base=2,padto=self._n))

		if sboxFt[0] != 0:
			sboxFt = [ g^^sboxFt[0] for g in sboxFt ]

		for i in xrange(self._n):
			x=1<<i
			M.set_column(i,ZZ(sboxFt[x]).digits(base=2,padto=self._n))

		sbox = range(self._length)

		V=ZZ(V.list(),2)

		for i in xrange(self._length):
			sbox[i]=vector(GF(2),ZZ(i).digits(base=2,padto=self._n))
			sbox[i]=M*sbox[i]
			sbox[i]=ZZ(sbox[i].list(),2)
			
			sbox[i]= sbox[i]^^V

			sbox[i]=sboxG[sbox[i]]

		V=vector(V.digits(2,padto=self._n))

		#print "~"*10
		#print "M:\n{0}".format(M)
		#print "V:\n{0}".format(V)
		#print "~"*10

		#print "sboxF  = {0}".format(sboxF)
		#print "sbox   = {0}".format(sbox)

		if sboxF == sbox:
			return [M,V]
		else:
			return None

	def check_pEA3(sboxF,sboxG):
		r'''
			Return function M and V if sboxF(x) = sboxG(x) + M * x + V
		'''
		M = matrix(GF(2),nrows=self._m,ncols=self._n)

		sboxFt=sboxF[:]
		sboxGt=sboxG[:]

		for i in xrange(self._length):
			sboxFt[i]=sboxGt[i]^^sboxFt[i]

		V = vector(GF(2),ZZ(sboxFt[0]).digits(base=2,padto=self._n))

		if sboxFt[0] != 0:
			sboxFt = [ g^^sboxFt[0] for g in sboxFt ]

		for i in xrange(self._n):
			x=1<<i
			M.set_column(i,ZZ(sboxFt[x]).digits(base=2,padto=self._n))

		sbox = range(self._length)

		V=ZZ(V.list(),2)

		for i in xrange(self._length):
			sbox[i]=vector(GF(2),ZZ(i).digits(base=2,padto=self._n))
			sbox[i]=M*sbox[i]
			sbox[i]=ZZ(sbox[i].list(),2)
			
			sbox[i]= sboxG[i]^^sbox[i]^^V

		V=vector(V.digits(2,padto=self._n))

		#print "~"*10
		#print "M:\n{0}".format(M)
		#print "V:\n{0}".format(V)
		#print "~"*10

		#print "sboxF  = {0}".format(sboxF)
		#print "sbox   = {0}".format(sbox)

		if sboxF == sbox:
			return [M,V]
		else:
			return None

	def check_pEA4(F,G):
		r'''
			Return function M1, M3 and V if F(x) = M1 * G(x) + M3 * x + V
			
			F(x) = F'(x) + L1(x) + V1
			G(x) = G'(x) + L2(x) + V2

			F'(x) + L1(x) + V1 = M1 * G'(x) + M1 * L2(x) + M1 * V2 + M3 * x + V

			F'(x) = M1 * G'(x) and L1(x) + V1 = M1 * L2(x) + M1 * V2 + M3 * x + V

			For known M1:
			L1(x) + V1 + M1 * L2(x) + M1 * V2 = M3 * x + V
		'''
		M1 = matrix(GF(2),nrows=self._m,ncols=self._n)
		M3 = matrix(GF(2),nrows=self._m,ncols=self._n)

		polF = F
		polG = G

		V1 = polF.constant_coefficient()
		V2 = polG.constant_coefficient()

		polF += V1
		polG += V2

		V1 = V1.integer_representation()
		V2 = V2.integer_representation()

		polFc=polF.coeffs()
		polFc += [self._P("0") for i in xrange(self._length-len(polFc))]
		polGc=polG.coeffs()
		polGc += [self._P("0") for i in xrange(self._length-len(polGc))]

		L1 = zero_vector(self._length).list()
		L2 = zero_vector(self._length).list()

		for i in xrange(self._n):
			if polFc[1<<i] != 0:
				L1[1<<i] = polFc[1<<i]
				polFc[1<<i] = 0

			if polGc[1<<i] != 0:
				L2[1<<i] = polGc[1<<i]
				polGc[1<<i]	= 0

		L1 = self._P(L1)
		L2 = self._P(L2)
		polF = self._P(polFc)
		polG = self._P(polGc)

		sboxF = range(self._length)
		sboxG = range(self._length)

		sboxL1 = [L1.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]
		sboxL2 = [L2.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]
		sboxF  = [polF.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]
		sboxG  = [polG.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]

		sboxFt=sboxF[:]
		sboxGt=sboxG[:]


		if len(set(sboxG).intersection(set([2^g for g in xrange(self._n)]))) != self._n:
			#print ">>> sboxG hasn't all values of {0} <<<".format([2^g for g in xrange(self._n)])
			return None

		# start find M1
		for i in xrange(self._n):
			x=sboxGt.index(1<<i)
			M1.set_column(i,ZZ(sboxFt[x]).digits(base=2,padto=self._n))

		sboxM = range(self._length)

		V = ZZ((M1*vector(GF(2),ZZ(V2).digits(base=2,padto=self._n))).list(),2) ^^ V1
		for i in xrange(self._length):
			sboxM[i] = sboxL1[i] ^^ ZZ((M1*vector(GF(2),ZZ(sboxL2[i]).digits(base=2,padto=self._n))).list(),2) ^^ V

		sboxT=sboxM[:]

		V = vector(GF(2),ZZ(sboxT[0]).digits(base=2,padto=self._n))
		
		if sboxT[0] != 0:
			sboxT = [ g^^sboxT[0] for g in sboxT ]

		# start find M3
		for i in xrange(self._n):
			x=1<<i
			M3.set_column(i,ZZ(sboxT[x]).digits(base=2,padto=self._n))

		sbox = range(self._length)

		sF  = [F.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]
		sG  = [G.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]		

		# start checking
		for i in xrange(self._length):
			sbox[i]=vector(GF(2),ZZ(sG[i]).digits(base=2,padto=self._n))

			sbox[i]=M1*sbox[i]

			tx=M3*vector(GF(2),ZZ(i).digits(base=2,padto=self._m))
			
			sbox[i]=vector(GF(2),[ZZ(sbox[i].get(j)) ^^ ZZ(tx.get(j)) ^^ ZZ(V.get(j)) for j in xrange(len(sbox[i]))])

			sbox[i]=ZZ(sbox[i].list(),2)

		if sbox == sF:
			return [M1,M3,V]
		else:
			return None

	if F is None:
		raise TypeError("Function 'F' must be presented")

	if G is None:
		raise TypeError("Function 'G' must be presented")

	degF = max([ZZ(i).popcount() for i in F.exponents()])
	degG = max([ZZ(i).popcount() for i in G.exponents()])

	if degG != degF:
		return False

	sboxF = range(self._length)
	sboxG = range(self._length)
	
	sboxF = [F.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]
	sboxG = [G.subs(self._K(ZZ(i).digits(2))).integer_representation() for i in xrange(self._length)]

	M1 = None
	V1 = None
	M2 = None
	V2 = None
	M3 = None

	for v2 in xrange(self._length):
		sbox=sboxF[:]
		sbox=[sbox[g^^v2] for g in xrange(self._length)]
		ret=check_pEA1(sbox,sboxG)

		if ret != None:
			M1=ret[0]
			V1=ret[1]
			V2=vector(GF(2),ZZ(v2).digits(base=2,padto=self._n))
			if functions == True:
				return [M1,V1,M2,V2,M3]
			else:
				return True

	for v1 in xrange(self._length):
		sbox=sboxF[:]
		sbox=[sbox[g]^^v1 for g in xrange(self._length)]
		ret=check_pEA2(sbox,sboxG)

		if ret != None:
			M2=ret[0]
			V2=ret[1]
			V1=vector(GF(2),ZZ(v1).digits(base=2,padto=self._n))
			if functions == True:
				return [M1,V1,M2,V2,M3]
			else:
				return True


	for v2 in xrange(self._length):
		sbox=sboxG[:]
		sbox=[sbox[g^^v2] for g in xrange(self._length)]
		ret=check_pEA3(sboxF,sbox)

		if ret != None:
			M3=ret[0]
			V1=ret[1]
			V2=vector(GF(2),ZZ(v2).digits(base=2,padto=self._n))
			if functions == True:
				return [M1,V1,M2,V2,M3]
			else:
				return True

	for v2 in xrange(self._length):
		polG=G.subs(self._P("x+{0}".format(self._K(ZZ(v2).digits(2))))).mod(self._P("x^{0}+x".format(self._length)))

		ret=check_pEA4(F,polG)

		if ret != None:
			M1=ret[0]
			M3=ret[1]
			V1=ret[2]
			V2=vector(GF(2),ZZ(v2).digits(base=2,padto=self._n))
			if functions == True:
				return [M1,V1,M2,V2,M3]
			else:
				return True

	return False


def cr_nonlinearity(self):
	r"""
		bits   time
		  10      7
		  11     25
		  12     97
		  13    405
		  14   1659
	"""
	return c_nonlinearity(self._S,1<<self._m)
