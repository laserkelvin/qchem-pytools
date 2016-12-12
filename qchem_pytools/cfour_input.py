
import os
import subprocess
import sys

class cfour_zmat:
	settings = {
		"basis": "ANO0",
		"method": "CCSD",
		"scf_conv": "9",
		"cc_conv": "9",
		"geo_conv": "5",
		"lineq_conv": "9",
		"zmat": """
		""",
		"frozen_core": "OFF",
		"comment": "Generated CFOUR ZMAT",
		"cfour_path": "/opt/cfour/cfour15/bin/",
	}

	zmat_template = """{comment}
{zmat}
*CFOUR(CALC_LEVEL={method}
BASIS={basis}
SCF_CONV={scf_conv}
CC_CONV={cc_conv}
GEO_CONV={geo_conv}
LINEQ_CONV={lineq_conv}
ABCDTYPE=AOBASIS
FROZEN_CORE={frozen_core}"""


	def __init__(self, comment, zmat, filename):
		self.settings["comment"] = comment
		self.settings["filename"] = filename
		self.settings["zmat"] = zmat
		self.settings["root_path"] = os.path.abspath(os.curdir)
		self.settings["file_path"] = self.settings["root_path"] + "/calcs/" + filename + "/"
		if self.settings["cfour_path"] not in os.environ["PATH"]:
			os.environ["PATH"] += ":" + self.settings["cfour_path"]
		self.output = self.settings["file_path"] + filename + ".out"


	def write_zmat(self):
		try:
			os.mkdir(self.settings["file_path"])
		except FileExistsError:
			cont = input("Folder already exists, continue? Y/N\t")
			if cont == "Y" or "y":
				with open(self.settings["file_path"] + "ZMAT", "w+") as WriteFile:
					self.input = self.input + ")\n"
					WriteFile.write(self.input)
			elif cont == "N" or "n":
				pass


	def build_zmat(self, extra=None):
		""" Method that will construct the ZMAT file. 
			Extra arguments should be provided as a list of strings;
			["VIB=EXACT", "EXCITE=EOMIP"]
		"""
		self.input = self.zmat_template.format_map(self.settings)
		if extra != None:
			for item in extra:
				self.input = self.input + item + "\n"


	def run_calc(self):
		command = self.settings["cfour_path"] + "xcfour"
		os.chdir(self.settings["file_path"])
		with open(self.output, "w+") as WriteFile:
			subprocess.run([command], stdout=WriteFile, stderr=WriteFile)
		os.chdir(self.settings["root_path"])
		print("Calculation finished.")


	def clean(self):
		os.chdir(self.settings["file_path"])
		subprocess.run(["xclean"])
		files = [
		"FCM*", "MOL", "OLDMOS", 
		"THETA", "OPTARC", "VPOUT", 
		"dens.dat", "basinfo.dat",
		"EFG",
		]
		for file in files:
			os.system("rm -r " + file)
		os.chdir(self.settings["root_path"])
