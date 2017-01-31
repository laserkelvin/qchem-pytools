
import os
import subprocess
import json
import time

class cfour_zmat:

    zmat_template = """{comment}
{zmat}
*CFOUR(CALC_LEVEL={method}
BASIS={basis}
SCF_CONV={scf_conv}
CC_PROGRAM={cc_program}
CC_CONV={cc_conv}
GEO_CONV={geo_conv}
LINEQ_CONV={lineq_conv}
MEM_UNIT=GB
MEMORY_SIZE={memory}
CHARGE={charge}
ABCDTYPE=AOBASIS
EXCITE={excite}
DBOC={dboc}
RELATIVISTIC={relativistic}
VIB={vib}
FREQ_ALGORITHM={freq_algorithm}
REFERENCE={reference}
MULTIPLICITY={multiplicity}
SCF_DAMPING={scf_damping}
FROZEN_CORE={frozen_core})
{footer}
"""

    def __init__(self, comment, zmat, filename, rerun=False):
        self.settings = {
            "basis": "ANO0",
            "method": "CCSD",
            "cc_program": "ECC",
            "memory": "1",
            "scf_conv": "9",
            "cc_conv": "9",
            "geo_conv": "5",
            "lineq_conv": "9",
            "reference": "RHF",
            "charge": "0",
            "multiplicity": "1",
            "excite": "NONE",
            "scf_damping": "1000",
            "dboc": "OFF",
            "relativistic": "OFF",
            "vib": "NO",
            "freq_algorithm": "ANALYTIC",
            "zmat": """
            """,
            "frozen_core": "OFF",
            "comment": "Generated CFOUR ZMAT",
            "cfour_path": "/opt/cfour/cfour15/bin/",
            "footer": "",
            "timestamp": time.strftime("%d/%m/%Y") + "\t" + time.strftime("%H:%M:%S")
        }
        self.rerun = rerun
        self.settings["comment"] = comment
        self.settings["filename"] = filename
        self.settings["zmat"] = zmat
        self.settings["root_path"] = os.path.abspath(os.curdir)
        self.settings["file_path"] = self.settings["root_path"] + "/calcs/" + filename + "/"
        if self.settings["cfour_path"] not in os.environ["PATH"]:
            os.environ["PATH"] += ":" + self.settings["cfour_path"]
        self.output = self.settings["file_path"] + filename + ".out"

    def write_zmat(self):
        if (os.path.exists(self.settings["file_path"])) is True:
            pass
        else:
            os.mkdir(self.settings["file_path"])
        if os.path.isfile(self.settings["file_path"] + "ZMAT"):
            cont = input("Folder already exists, continue? Y/N\t")
            if cont == "Y" or "y":
                with open(self.settings["file_path"] + "ZMAT", "w+") as WriteFile:
                    WriteFile.write(self.input)
            elif cont == "N" or "n":
                pass
        else:
            with open(self.settings["file_path"] + "ZMAT", "w+") as WriteFile:
                WriteFile.write(self.input)

    def build_zmat(self, extra=None):
        """ Method that will construct the ZMAT file.
            Extra arguments should be provided as a list of strings;
            ["VIB=EXACT", "EXCITE=EOMIP"]
        """
        self.input = self.zmat_template.format_map(self.settings)
        if extra is not None:
            for item in extra:
                self.input = self.input + item + "\n"
        print(self.input)

    def run_calc(self):
        command = self.settings["cfour_path"] + "xcfour"
        os.chdir(self.settings["file_path"])
        if self.rerun is True:
            with open(self.output, "w+") as WriteFile:
                subprocess.run([command], stdout=WriteFile, stderr=WriteFile)
                print("Calculation finished.")
                self.clean()
        elif self.rerun is False:
            print("Rerun is set to false, calculation not run.")
        os.chdir(self.settings["root_path"])

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

    def save_json(self, Filename=None):
        if Filename is None:
            Filename = "json/" + self.settings["filename"] + ".settings.json"
        with open(Filename, "w+") as WriteFile:
            json.dump(self.settings, WriteFile)
