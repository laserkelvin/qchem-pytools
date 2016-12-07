""" This is a crudely written class for parsing
    output from the electronic structure package "CFOUR".

    The idea behind this set of routines is to be able to
    quickly extract all of the useful information out of
    a CFOUR calculation. The information is stored as a Python
    dictionary, which has the syntax of JSON i.e. portable.

    The idea is to also be able to generate quick HTML
    reports so that I can file them into Jupyter notebooks.

    Some of the code is written pythonically, while
    others are not so much...

"""

class OutputFile:
    InfoDict = {
        "filename": " ",
        "basis": " ",
        "success": False,
        "method": " ",
        "dipole moment": [],
        "point group": " ",
        "energies": {
            "final_energy": 0.,
            "scf_cycles": dict(),
            "cc_cycles": dict(),
        },
        "coordinates": [],
        "input zmat": [],
        "final zmat": [],
        "natoms": 0,
        "gradient norm": dict(),
    }
    def __init__(self, File):
        self.InfoDict["filename"] = File
        self.parse()


    def parse(self):
        scf_iter = 0
        scf_cycle = 0
        cc_counter = 0
        geo_counter = 0
        skip_counter = 0
        CurrentCoords = []
        DipoleFlag = False
        RotFlag = False
        SCFFlag = False
        CCFlag = False
        IZMATFlag = False        # Initial ZMAT file
        FZMATFlag = False        # Final ZMAT file
        ReadCoords = False
        with open(self.InfoDict["filename"], "r") as ReadFile:
            for LineIndex, Line in enumerate(ReadFile):
                if ("The final electronic energy is") in Line:
                    ReadLine = Line.split()
                    self.InfoDict["energies"]["final_energy"] = float(ReadLine[5])
                    self.InfoDict["success"] = True
                if ("The full molecular point group ") in Line:
                    ReadLine = Line.split()
                    self.InfoDict["point group"] = ReadLine[6]
                if ("BASIS=") in Line:
                    ReadLine = Line.split("=")
                    self.InfoDict["basis"] = ReadLine[1].split()[0]
                if ("CALC_LEVEL") in Line:
                    ReadLine = Line.split("=")
                    self.InfoDict["method"] = ReadLine[1].split()[0]
                if ("EXCITE=") in Line:
                    ReadLine = Line.split("=")
                    self.InfoDict["method"] += "-" + ReadLine[1].split()[0]
                if RotFlag is True:            # if flagged to read the rotational constants
                    ReadLine = Line.split()
                    self.InfoDict["rotational constants"] = [float(Value) for Value in ReadLine]
                    RotFlag = False
                if ("Rotational constants (in MHz)") in Line:
                    RotFlag = True
                if SCFFlag is True:
                    ReadLine = Line.split()
                    if len(ReadLine) == 3:
                        ReadLine = [Item.replace("D", "E") for Item in ReadLine]
                        try:
                            CurrentSCF[scf_iter] = [float(Value) for Value in ReadLine]
                            scf_iter += 1
                        except ValueError:
                            pass
                if ("Total Energy") in Line:
                    SCFFlag = True
                    scf_iter = 0
                    CurrentSCF = dict()
                if ("SCF has converged") in Line:
                    self.InfoDict["energies"]["scf_cycles"][scf_cycle] = CurrentSCF
                    SCFFlag = False
                    scf_cycle += 1
                if FZMATFlag is True:
                    if ("********") in Line:
                        skip_counter += 1
                    elif skip_counter == 1:
                        self.InfoDict["final zmat"].append(Line)
                    elif skip_counter == 2:
                        FZMATFlag = False
                        skip_counter = 0
                if ("Final ZMATnew file") in Line:
                    FZMATFlag = True
                if ReadCoords is True:
                    if ("----------") in Line:
                        skip_counter += 1
                    elif skip_counter == 1:
                        CurrentCoords.append(Line.split())
                    elif skip_counter == 2:
                        self.InfoDict["coordinates"] = CurrentCoords
                        ReadCoords = False
                        CurrentCoords = []
                        skip_counter = 0
                if ("Coordinates (in bohr)") in Line:
                    ReadCoords = True
                if ("Conversion factor used") in Line:
                    DipoleFlag = False
                if DipoleFlag is True:
                    ReadLine = Line.split()
                    if len(ReadLine) == 3:
                        Dipole.append(ReadLine)
                if ("au             Debye") in Line:
                    Dipole = []
                    DipoleFlag = True
                if ("Molecular gradient norm") in Line:
                    ReadLine = Line.split()
                    self.InfoDict["gradient norm"][geo_counter] = float(ReadLine[3])
                    geo_counter += 1
                if CCFlag is True:
                    if skip_counter == 2:
                        self.InfoDict["energies"]["cc_cycles"][cc_counter] = CCCycle
                        CCFlag = False
                        cc_counter += 1
                    elif ("-------") in Line:
                        skip_counter += 1
                    else:
                        ReadLine = Line.split()[:3]
                        CCCycle.append([float(Value) for Value in ReadLine])
                if ("Iteration         Energy              Energy") in Line:
                    skip_counter = 0
                    CCCycle = []
                    CCFlag = True
        self.InfoDict["natoms"] = len(self.InfoDict["coordinates"])


    def print_results(self):
        return self.InfoDict


    def export_xyz(self, Filename):
        with open(Filename, "w+") as WriteFile:
            WriteFile.write("Output geometry for " + self.InfoDict["filename"] + " in Bohr\n")
            WriteFile.write(str(self.InfoDict["natoms"]))
            WriteFile.write("\n")
            for Line in self.InfoDict["coordinates"]:
                for Piece in Line:
                    WriteFile.write(Piece)
                    WriteFile.write(" ")
                WriteFile.write("\n")


    def print_xyz(self):
        for Line in self.InfoDict["coordinates"]:
            print(Line)


    #def save_html(self, Filename=None):
        #T = Template(HTMLTemplate)
        #C = Context(self.InfoDict)
        #HTMLOut = T.render(C)
        #if Filename is None:
        #    Filename = self.InfoDict["filename"] + "_report.html"
        #with open(Filename, "w+") as WriteFile:
        #    WriteFile.writelines(HTMLOut)


    def save_json(self, Filename=None):
        import json
        if Filename is None:
            Filename = self.InfoDict["filename"] + ".json"
        with open(Filename, "w+") as WriteFile:
            json.dump(self.InfoDict, WriteFile)

