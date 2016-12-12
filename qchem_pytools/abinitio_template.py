""" Generate a template for an electronic structure program
	input file.

"""

def generate_abinitio_input(Package="CFOUR"):
	input = ""
	if Package == "CFOUR":
		input = '''
		*CFOUR(CALC_LEVEL=
		BASIS=
		FROZEN_CORE=
		CC_CONV=9
		SCF_CONV=9
		LINEQ_CONV=9
		ESTATE_CONV=9
		GEO_CONV=5
		ABCDTYPE=AOBASIS
		)
		'''
	return input