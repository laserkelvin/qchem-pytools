""" In the spirit of reproducibility, these routines will automatically
	generate a file directory tree for performing calculations on a given
	molecule.

	The idea will be to have one notebook per "experiments", which will
	collate all of the results and produce some interpretation based on
	whatever has been done.

"""

import os

top_level = [
	"figures",        # from the notebook outputs
	"suppinfo",       # whatever extra data required
	"calcs",          # store the calculation data
	"json",        # collated information from calculations
	"docs"            # possibly reports of the analysis
]

navigation = {
	"top": os.path.abspath("./")
}


def setup_folders():
	for folder in top_level:
		try:
			os.mkdir(folder)
		except FileExistsError:
			print(folder + " already exists.")
			pass
