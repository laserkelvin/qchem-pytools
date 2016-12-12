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

import pandas as pd
import numpy as np
import os
from qchem_pytools import figure_settings
from qchem_pytools import html_template

class OutputFile:
    InfoDict = {
        "filename": " ",
        "basis": " ",
        "success": False,
        "method": " ",
        "dipole moment": [0., 0., 0.],
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
        "nscf": 0.,
        "ncc": 0.,
        "avg_scf": 0.,
        "avg_cc": 0.,
        "gradient norm": [],
        "paths": {
           "root": " ",
            "json": " ",
            "figures": " ",
            "calcs": " ",
            "output": " ",
        }
    }


    def __init__(self, File):
        for key in self.InfoDict["paths"]:
            self.InfoDict["paths"][key] = os.path.abspath("./" + key) + "/"
        self.InfoDict["paths"]["full"] = self.InfoDict["filename"] = File
        self.InfoDict["filename"] = os.path.split(File)[1].split(".")[0]
        self.parse()
        self.plot_generation()
        self.save_json(self.InfoDict["paths"]["json"] + self.InfoDict["filename"] + ".json")
        html_template.html_report(self.InfoDict)


    def parse(self):
        scf_iter = 0
        scf_cycle = 0
        cc_cycle = 0
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
        with open(self.InfoDict["paths"]["full"], "r") as ReadFile:
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
                            CurrentSCF.append(float(ReadLine[1]))      # Take the energy
                            scf_iter += 1
                        except ValueError:
                            pass
                if ("Total Energy") in Line:
                    SCFFlag = True
                    scf_iter = 0
                    CurrentSCF = []
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
                        ReadLine = Line.split()
                        CurrentCoords.append([ReadLine[0], ReadLine[2], ReadLine[3], ReadLine[4]])
                    elif skip_counter == 2:
                        self.InfoDict["coordinates"] = CurrentCoords
                        ReadCoords = False
                        CurrentCoords = []
                        skip_counter = 0
                if ("Coordinates (in bohr)") in Line:
                    skip_counter = 0
                    ReadCoords = True
                if ("Conversion factor used") in Line:
                    self.InfoDict["dipole moment"] = Dipole
                    DipoleFlag = False
                if DipoleFlag is True:
                    ReadLine = Line.split()
                    if len(ReadLine) == 3:
                        Index = ["x", "y", "z"].index(ReadLine[0])
                        Dipole[Index] = float(ReadLine[2])
                if ("au             Debye") in Line:
                    Dipole = [0., 0., 0.]
                    DipoleFlag = True
                if ("Molecular gradient norm") in Line:
                    ReadLine = Line.split()
                    self.InfoDict["gradient norm"].append(float(ReadLine[3]))
                    geo_counter += 1
                if CCFlag is True:
                    if skip_counter == 2:
                        self.InfoDict["energies"]["cc_cycles"][cc_cycle] = CurrentCC
                        CCFlag = False
                        cc_cycle += 1
                    elif ("-------") in Line:
                        skip_counter += 1
                    else:
                        ReadLine = Line.split()[:3]
                        CurrentCC.append([float(ReadLine[1]), float(ReadLine[2])])
                if ("Iteration         Energy              Energy") in Line:
                    skip_counter = 0
                    CurrentCC = []
                    CCFlag = True
        self.InfoDict["natoms"] = len(self.InfoDict["coordinates"])
        self.InfoDict["nscf"] = len(self.InfoDict["energies"]["scf_cycles"])
        self.InfoDict["ncc"] = len(self.InfoDict["energies"]["cc_cycles"])

        scf_iteration_data = [len(self.InfoDict["energies"]["scf_cycles"][cycle]) for cycle in self.InfoDict["energies"]["scf_cycles"]]
        cc_iteration_data = [len(self.InfoDict["energies"]["cc_cycles"][cycle]) for cycle in self.InfoDict["energies"]["cc_cycles"]]
        self.InfoDict["avg_scf"] = np.average(scf_iteration_data)
        self.InfoDict["avg_cc"] = np.average(cc_iteration_data)


    def print_results(self):
        return self.InfoDict


    def plot_generation(self):
        """ Method that will generate a summary of energy convergence """
        first_scf = pd.DataFrame(
                data=self.InfoDict["energies"]["scf_cycles"][0],
                columns=["SCF energy"]
            )
        last_scf = pd.DataFrame(
                data=self.InfoDict["energies"]["scf_cycles"][self.InfoDict["nscf"] - 1],
                columns=["SCF energy"]
            )
        scf_figure = figure_settings.GenerateSubPlotObject(
                ["First SCF cycle", "Final SCF cycle"],
                1,
                2,
                "Landscape",
            )
        for index, plot in enumerate([first_scf, last_scf]):
            plot_object = figure_settings.DefaultScatterObject(
                    plot.index,
                    plot["SCF energy"]
                )
            scf_figure.append_trace(plot_object, 1, index + 1)
            scf_figure["layout"]["xaxis" + str(index + 1)].update(
                    {
                        "title": "Iterations"
                    }
                )
        scf_figure["layout"].update(
                width= 900.,
                height= (900. / 1.6),
                showlegend=False
            )
        scf_figure["layout"]["yaxis"].update(
                {
                    "title": "Energy (Ha)"
                }
                )
        with open(self.InfoDict["paths"]["figures"] + self.InfoDict["filename"] + ".scf_report.html", "w+") as WriteFile:
            WriteFile.write(
                figure_settings.plot(
                        scf_figure,
                        output_type="div",
                        auto_open=False
                    )
            )
        if len(self.InfoDict["gradient norm"]) > 0:
            geometry_opt = pd.DataFrame(
                    data=self.InfoDict["gradient norm"]
                )
            geometry_plot = [
                figure_settings.DefaultScatterObject(
                        X=geometry_opt.index,
                        Y=geometry_opt[0],
                        Name="Molecular gradient",
                    )
            ]
            geometry_layout = figure_settings.DefaultLayoutSettings()
            geometry_layout["title"] = "Geometry optimisation progress"
            geometry_layout["xaxis"]["title"] = "Iteration"
            geometry_layout["yaxis"]["title"] = "Molecular gradient norm"
            with open(self.InfoDict["paths"]["figures"] + self.InfoDict["filename"] + ".geo_report.html", "w+") as WriteFile:
                WriteFile.write(
                        figure_settings.plot(
                                {"data": geometry_plot, "layout": geometry_layout},
                                output_type="div",
                                auto_open=False
                            )
                    )


    def export_xyz(self, Filename=None):
        if Filename == None:
            Filename = self.InfoDict["filename"] + ".xyz"
        Filename = self.path["root"] + "/" + Filename
        with open(Filename, "w+") as WriteFile:
            WriteFile.write(str(self.InfoDict["natoms"]))
            WriteFile.write("Output geometry for " + self.InfoDict["filename"] + " in Bohr\n")
            WriteFile.write("\n")
            for Line in self.InfoDict["coordinates"]:
                for Piece in Line:
                    WriteFile.write(Piece)
                    WriteFile.write(" ")
                WriteFile.write("\n")


    def print_xyz(self):
        for Line in self.InfoDict["coordinates"]:
            print(Line)


    def save_json(self, Filename=None):
        import json
        if Filename is None:
            Filename = self.InfoDict["filename"] + ".json"
        with open(Filename, "w+") as WriteFile:
            json.dump(self.InfoDict, WriteFile)
