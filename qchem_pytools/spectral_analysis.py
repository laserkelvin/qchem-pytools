#!/bin/python

from collections import OrderedDict
import inspect
import numpy as np
import pandas as pd
from qchem_pytools import notebook_tools as NT
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import constants
import os
import peakutils
from scipy import fftpack
from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit, leastsq, least_squares
from bokeh.palettes import brewer
from bokeh.plotting import figure, show
from numba import jit
import seaborn as sns

# SpectralAnalysis.py

# Module containing functions for day-to-day spectral analysis, including re-packaged
# routines from FittingRoutines.py
# My plan is to make FittingRoutines obsolete, as it has A LOT of poorly written
# routines.

###################################################################################################

""" General Plot settings """

matplotlib.style.use(['seaborn-pastel', 'fivethirtyeight'])      # Theme
D3State = False

###################################################################################################

""" Constants and small persistent dictionaries """

""" Values for oxygen calibration as a dictionary
OxygenKER is the total translational energy release for the molecule,
while OxygenAtomSpeed is self explanatory.
 """
OxygenKER = OrderedDict({"3P": 0.21585*2,
                         "5P": 0.33985*2,
                         "3S": 0.95*2,
                         "5S": 1.14*2}
                       )

OxygenAtomSpeed = OrderedDict({"3P": 1141.016103985279,
                               "5P": 1431.7242621001822,
                               "3S": 2393.7431587974975,
                               "5S": 2622.2142498941203}
                             )

# Unit conversion, kB from J/mol to 1/(cm mol)
kcm = constants.physical_constants["Boltzmann constant in inverse meters per kelvin"][0] / 100.0

###################################################################################################

""" Classes """

class Spectrum:
    instances = []
    def __init__(self, File=None, CalConstant=1., Reference=None):
        """ If nothing is supplied, we'll just initialise a new
        instance of Spectrum without giving it any data.
        That way I don't HAVE to load a file.
        """
        self.CalibrationConstant = CalConstant
        self.CalibratedWavelengths = [0., 0.,]
        self.PumpWavelength = 0.
        self.ProbeWavelength = 0.
        self.Annotations = dict()
        self.Comments = ""
        if Reference is not None:                       # If we give it a logbook reference
            self.Reference = Reference
        self.Labels = dict()
        if File is not None:
            self.Data = LoadSpectrum(File, CalConstant)
            self.PlotLabels()                       # Only set up plotlabels when we have data

        Spectrum.instances.append(self)

    def CalibrateWavelengths(self, Wavelengths):
        """ Wavelengths given in as 2-tuple list
        and sets the Dataframe index to the calibrated
        wavelengths
        """
        NData = len(self.Data.index)
        NewAxis = np.linspace(num=NData, *Wavelengths)
        self.Data.index = NewAxis

    ###################################### Spectrum I/O ######################################

    def AddData(self, NewData, Name):
        self.Data[Name] = NewData

    def DeleteData(self, Name):
        del self.Data[Name]

    def ExportData(self, CustomSuffix=False, Suffix=None):
        try:
            os.mkdir("DataExport")
        except OSError:
            pass
        if CustomSuffix == False:
            try:
                FilePath = "./DataExport/" + self.Reference + "_export.csv"
            except AttributeError:
                FilePath = raw_input(" No reference found, please specify file.")
        elif CustomSuffix == True and Suffix == None:
            try:
                Suffix = raw_input("Please enter a suffix to be used for export, e.g. _export")
                FilePath = "./DataExport/" + self.Reference + Suffix
            except AttributeError:
                Reference = raw_input("No reference found, please provide one.")
                FilePath = "./DataExport/" + Reference + Suffix
        elif Suffix != None:
            try:
                FilePath = self.Reference + Suffix
            except AttributeError:
                Reference = raw_input("No reference found, please provide one.")
                FilePath = "./DataExport/" + Reference + Suffix
        self.Data.to_csv(FilePath)
        print(" File saved to:\t" + FilePath)

    def KERfromPES(self, File, Mass, Units="cm"):
        self.PES = LoadSpectrum(File, self.CalibrationConstant)
        self.Data = ConvertSpeedDistribution(self.PES, Mass, Units=Units)

    def ExportFits(self, Suffix="_fit.csv"):
        try:
            os.mkdir("FittingResults")
        except OSError:
            pass
        try:
            self.FitResults.to_csv("./FittingResults/" + self.Reference + Suffix)
        except AttributeError:
            Reference = raw_input(" No reference found, give me a name.")
            self.FitResults.to_csv("./FittingResults/" + Reference + Suffix)
            print(" File saved to:\t" + Reference + Suffix)

    def LoadDatabaseSpectrum(self, Database, Reference, Bogscan=False):
        Data = NT.LoadReference(Database, Reference, Verbose=False)
        Options = dict()
        if len(Data.keys()) != 1:       # Only ask if there's more than one file
            for Index, Key in enumerate(Data.keys()):
                Options[Index] = Key
            print(Options)
            Selector = int(raw_input(" Please specify which key to load"))
            Filename = Options[Selector]
        else:
            Filename = Data.keys()[0]
        self.Data = pd.DataFrame(data=Data[Filename])
        if Bogscan is True:
            if np.float(self.Data[1].sum()) > 0 is not True:    # Why does this not work
                self.Data.index = np.array(self.Data[1])
            else:
                self.Data.index = np.array(self.Data[0])
            del self.Data[0]
            del self.Data[1]
            del self.Data[3]
            self.Data.columns = ["Y Range"]
            if self.Data["Y Range"].sum() > -10. is True:
                self.Data["Y Range"] = self.Data["Y Range"]     # Make intensity positive

    def ReadBogScan(self, File):
        """ Special function for reading data files from
        BogScan
        """
        DataFrame = pd.read_csv(File, delimiter="\t", header=None)
        if DataFrame[1].sum > 0. == True:                   # If we've got calibrated wavelengths
            X = np.array(DataFrame[1])
            print(" Using calibrated wavelengths")
        else:
            X = np.array(DataFrame[0])
            print(" Using bogscan wavelengths")
        if DataFrame[2].sum > -100. == True:                # if the axis is negative
            Y = np.array(DataFrame[2])                      # I like plotting positive!
        else:
            Y = -np.array(DataFrame[2])
        NewDataFrame = FormatData(X, Y)
        self.Data = NewDataFrame

    ###################################### Plotting ######################################

    def Plot(self, Column="Y Range", Labels=None, Legend=True, Interface="pyplot"):
        self.PlotLabels()
        if Labels is not None:
            self.SetLabels(Labels)
        Target = FormatData(self.Data.index, np.array(self.Data[Column]))
        PlotData(DataFrame=Target, Labels=self.Labels, Interface=Interface, Legend=Legend)

    def PlotAll(self, Columns=None, Labels=None, Legend=True, Interface="pyplot", PlotTypes=None, ShowAnnotations=True):
        self.PlotLabels()
        Annotations = None
        if Labels is not None:
            self.SetLabels(Labels)
        if PlotTypes is not None:
            self.PlotTypes = PlotTypes
        if ShowAnnotations is True:
            Annotations = self.Annotations
        PlotData(DataFrame=self.Data,
                 Columns=Columns,
                 Labels=self.Labels,
                 Interface=Interface,
                 Legend=Legend,
                 PlotTypes=PlotTypes,
                 Annotations=Annotations
                 )

    def PlotLabels(self, Column=None):
        if Column is None:
            Column = self.Data.keys()[0]       # if no data is used to initialise
        Labels = {"X Label": "X Axis",
                      "Y Label": "Y Axis",
                      "Title": self.Reference,
                      "X Limits": [min(self.Data.index), max(self.Data.index)],
                      "Y Limits": [min(self.Data[Column]), max(self.Data[Column]) + max(self.Data[Column]) * 0.1],
                     }
        self.SetLabels(Labels)

    def SetLabels(self, Labels):
        """ For setting the labels for a plot
            If done this way, we can keep whatever is already set
            unless changed.
        """
        for Key in Labels:
            self.Labels[Key] = Labels[Key]

    ###################################### Analysis ######################################

    def Fit(self, Model, Column="Y Range", Interface="pyplot", Verbose=True):
        """ Calls the FitModel function to fit the Data contained in this
        instance.

        Requires Model instance reference as input, and will only fit the
        column labelled "Y Range" in data (i.e. the initially loaded data)

        Sets attributes of instance corresponding to the optimised parameters,
        the fit report, the fitted curves dataframe and covariance matrix

        Verbose mode will toggle between lots of printing or no printing at all
        when you may want to do batch fits (like in BootstrapAnalysis).
        """
        try:
            FittingData = FormatData(self.Data.index, np.array(self.Data[Column]))
            self.Opt, self.Report, self.FitResults, self.Cov = FitModel(FittingData, Model, Verbose=Verbose)
            #self.Report["Errors"] = self.Cov.diagonal()
            #self.Report.columns = ["Values", "Errors"]       # pack parameters and errors
            if Verbose == True:
                PlotData(self.FitResults, Labels=self.Labels, Interface=Interface)
            elif Verbose == False:
                pass
            self.Data["Fit Results"] = self.FitResults["Model Fit"]
            self.Data["Residuals"] = self.FitResults["Residuals"]
        except RuntimeError:
            print(''' No data column eligible for fitting/"Y Range" doesn't exist. ''')
            print(''' You may need to repack the data manually using FormatData. ''')

    def DetectPeaks(self, Column=None, Threshold=0.3, MinimumDistance=30.):
        """ Calls the peak finding function from peakutils that will
        sniff up peaks in a spectrum. This class method will then
        store that information as an attribute
        """
        PeakIndices = PeakFinding(self.Data, Column, Threshold, MinimumDistance)
        PeakX = [self.Data.index[Index] for Index in PeakIndices]
        PeakY = [1.2 for Peak in PeakIndices]
        self.Peaks = {"Peak Indices": PeakIndices,
                      "Threshold": Threshold,
                      "Minimum Distance": MinimumDistance}
        for PeakNumber, PeakIndex in enumerate(PeakIndices):
            self.Annotations[PeakNumber] = {"type": "vline",
                                            "position": self.Data.index[PeakIndex],
                                            "text": str(PeakNumber)}
        try:
            del self.PeakReport
        except AttributeError:
            pass
        self.PeakReport = pd.DataFrame(data=list(zip(PeakIndices, PeakX, PeakY)), columns=["Indices", "X Value", "Stick"])

    def Smooth(self, WindowSize=5, Column="Y Range"):
        """ Smooths the experimental data using a Svatizky-Golay filter
            with the specified window size.
        """
        try:
            self.Data[Column + "-Smoothed"] = np.array(NT.SGFilter(self.Data[Column], WindowSize))
        except KeyError:
            print(" No column " + Column + " found in DataFrame.")

    def NormaliseColumn(self, Columns=["Y Range"]):
        for Column in Columns:
            if "-Normalised" not in Column:
                NormaliseColumn(self.Data, Column=Column)
            else:
                pass

    def CalcAvailableEnergy(self, DissociationEnergy):
        """ Dissociation energy is in 1/cm, this routine calculates
            the available energy for a spectrum object
        """
        if self.PumpWavelength != 0.:
            self.Eavail = 1e7 / self.PumpWavelength - DissociationEnergy
        else:
            self.PumpWavelength = float(raw_input("No pump wavelength found, please specify one."))
            self.Eavail = 1e7 / self.PumpWavelength - DissociationEnergy
        return self.Eavail

    def BootstrapAnalysis(self, Model, DampFactor=0.5, Trials=100):
        """ Routine that will take a Spectral object and
            perform Bootstrap analysis: will generate synthetic data
            from a fit and repeatedly refit to get a sense of
            uncertainty
        """
        try:
            #OriginalFit = np.array(self.FitResults["Model Fit"])
            OriginalData = np.array(self.FitResults["Data"])
            X = np.array(self.FitResults.index)
        except AttributeError:
            print(" No fit results detected.")
            exit
        BootstrapResults = []
        for Trial in xrange(Trials):
            if ((float(Trial) / float(Trials)) * 100) % 20 == 0:      # check percentage progress
                print("Trials Done:\t" + str(Trial))
            NewData = AddNoise(OriginalData, DampFactor)
            try:
                NewFrame = FormatData(X, NewData)
                with NT.suppress_stdout():
                    Opt, Report, Curves, Cov = FitModel(NewFrame,
                                                        Model)
                BootstrapResults.append(Opt)
            except RuntimeError:
                pass
        self.BootstrapResults = pd.DataFrame(data=BootstrapResults,
                                            columns=Model.Variables)
        self.BootstrapReport = self.BootstrapResults.describe()
        self.BootstrapReport
        print(" Finished analysis. Check object BootstrapReport")

    def IntegrateColumn(self, Column="Y Range", Range=None):
        """ Return the integral of a column
            The range can be specified by a 2-tuple list by
            supplying the indices to begin and end.
        """
        X = np.nan_to_num(self.Data.index)              # Ensure there are no NaNs
        Y = np.nan_to_num(np.array(self.Data[Column]))  # in our arrays
        if Range is not None:
            X = X[Range[0]:Range[1]]
            Y = Y[Range[0]:Range[1]]
        return np.trapz(Y, X)

    def CalculateBranching(self):
        """ Calculate the branching fraction between dataframe items
            if Y Range is present, we assume it's the experimental
            data and treat that as the cumulative sum.
        """
        Integrals = dict()
        CumSum = 0.
        for Key in self.Data:
            Integrals[Key] = [self.IntegrateColumn(Column=Key)]
            CumSum = CumSum + self.IntegrateColumn(Column=Key)
        if "Y Range" in self.Data.keys():
            CumSum = self.IntegrateColumn(Column="Y Range")
        for Key in self.Data:
            Integrals[Key].append(Integrals[Key][0] / CumSum)
        self.Fractions = pd.DataFrame.from_dict(Integrals)
        self.Fractions.index = ["Integral", "Fraction"]
        return self.Fractions

###################################################################################################

class Model:
    """ Class for fitting models and functions. Ideally, I'll have a
    pickle file containing a lot of these model instances so that
    they'll be out of the way

    Input:
    Variables - Dictionary containing all of the variables required,
    with the keys named exactly as the function requests.

    ObjectiveFunction - Reference to a model function
    """
    def __init__(self, FunctionName):
        self.FunctionName = FunctionName
        self.Variables = {}
    def SetFunction(self, ObjectiveFunction, Verbose=True):
        """ Sets what the model functionw will be. Will also automagically
        initialise a dictionary that will hold all of the variables required,
        as well as the boundary conditions
        """
        self.Function = ObjectiveFunction
        self.Variables = OrderedDict.fromkeys(inspect.getargspec(ObjectiveFunction)[0])
        try:
            del self.Variables["x"]                  # X keeps getting picked up, get rid of it
            del self.Variables["X"]
        except KeyError:
            pass
        self.BoundaryConditions = ([-np.inf for Variable in self.Variables],
                                   [np.inf for Varaible in self.Variables])
        if Verbose == True:
            print(" Initialised variable dictionary:")
            print(self.Variables)
    def NewSetFunction(self, FunctionList, Verbose=True):
        self.FunctionList = FunctionList
        self.VariableDict = OrderedDict()
        for Function in FunctionList:
            self.VariableDict[Function.func_name] = OrderedDict.fromkeys(inspect.getargspec(Function)[0])
            try:
                del self.VariableDict[Function.func_name]["x"]
                del self.VariableDict[Function.func_name]["X"]
            except KeyError:
                pass
    def SetVariables(self, Variables, Verbose=True):
        self.Variables = UpdateDictionary(self.Variables, Variables)
        if Verbose == True:
            print(" Variables set to:")
            print(str(self.Variables))
    def SetBounds(self, Bounds, Verbose=True):
        """ Method for setting the boundary conditions for curve_fit,
        requires input as 2-tuple list ([], [])
        """
        self.BoundaryConditions = Bounds
        if Verbose == True:
            print(" Boundary conditions set to:" )
            print(str(self.BoundaryConditions))
    def ResetAttributes(self):
        """ Wipes the variable and boundary conditions
        """
        self.SetFunction(self.Function)

###################################################################################################

""" Data formatting and comprehension """

def ConvertSpeedDistribution(DataFrame, Mass, Units="cm"):
    """ Function to convert a data frame of speed distribution
    into a kinetic energy data frame. Requires the mass in amu
    as input, and returns the kinetic energy dataframe with the
    index as the KER, and "Y Range" as P(E).
    """
    Conversion = {"cm": 83.59,           # Dictionary of unit conversion from
                  "kJ": 1.,              # kJ/mol to 1/cm and eV
                  "eV": 0.0103636}
    KER = (np.square(DataFrame.index) * (Mass / 1000.) / 2000.) * Conversion[Units]
    PE = DataFrame["Y Range"].values / (Mass / 1000. * DataFrame.index)
    PE = [0. if Value < 0. else Value for Value in PE ]           # set negative values to zero
    KERDataFrame = pd.DataFrame(data=PE, index=KER, columns=["Y Range"])
    KERDataFrame = KERDataFrame.dropna(axis=0)
    #for index, value in enumerate(KERDataFrame)
    return KERDataFrame

def ConvertKERDistribution(DataFrame, Mass):
    Column = DataFrame.keys()[0]              # take the first Y set
    Speed = np.sqrt(((np.array(DataFrame.index) / 83.59) * 2000.) / (Mass / 1000.))
    PS = np.array(DataFrame[Column] * (Mass / 1000.)) * DataFrame.index
    SpeedDataFrame = pd.DataFrame(data=PS, index=Speed, columns=["Y Range"])
    SpeedDataFrame = SpeedDataFrame.dropna(axis=0)
    return SpeedDataFrame

def InvertInternalEnergy(DataFrame, MassFraction, Eavail):
    NewDF = DataFrame.copy()                                 # Copy the dataframe over
    NewDF.index = NewDF.index / MassFraction                 # Convert to TKER
    NewDF.index = Eavail - NewDF.index                       # Convert to Internal Energy
    return NewDF

def NormaliseColumn(DataFrame, Column="Y Range"):
    """ Routine to normalise a column in pandas dataframes
    """
    DataFrame[Column + "-Normalised"] = DataFrame[Column] / np.max(DataFrame[Column])

def Dict2List(Dictionary):
    List = [Dictionary[Item] for Item in Dictionary]
    return List

def GetFitParameters(FitReport, Parameters):
    """ Takes a input list of Parameter names
        and will return the optimised parameters
        in a dictionary
    """
    ParameterOut = dict()
    for Parameter in Parameters:
        ParameterOut[Parameter] = FitReport["Values"][Parameter]
    return ParameterOut

def UnpackDict(**args):
    """ I don't know if there's a better way of doing this,
    but I wrote this to upack a dictionary so we can parse
    a class dictionary and unpack it into curve_fit
    """
    print(args)

def ExtractFunctionParameters(Parameters, Report):
    """ Function for getting the parameters from a fit report
        produced by my fitting function.

        Takes a list of
        strings with names of the variables you want, and
        extracts them from the report
    """
    VariableDictionary = OrderedDict()
    for Variable in Parameters:
        VariableDictionary[Variable] = Report["Values"][Variable]
    return VariableDictionary

def FormatData(X, Y):
    """ Function to format data into a pandas data frame for
    fitting. In case I'm too lazy to set it up myself.
    """
    return pd.DataFrame(data=Y, columns=["Y Range"], index=X)

def DistributionStatistics(DataFrame):
    """ Calculate the numerical average for each column
        of a dataframe
    """
    Dictionary = {}
    X = np.array(DataFrame.index)
    for Key in DataFrame.keys():
        Y = np.array(DataFrame[Key])
        Y = Y / np.trapz(Y, X)                      # normalise to integral
        Expec = np.trapz(np.multiply(Y, X), X)      # expec = f(x) x dx
        Squares = np.multiply((X - Expec)**2, Y)
        StdDev = np.sqrt(np.trapz(Squares, X) / np.trapz(Y, X))
        Dictionary[Key] = {"Expec": Expec,
                           "StdDev": StdDev}
    return Dictionary

def SubtractSpectra(A, B):
    """ Takes input as two instances of Spectrum class, and does
    a nifty subtraction of the spectra A - B by interpolating
    B into the X axis range of A
    """
    XA = A.Data.index
    YA = A.Data.as_matrix(["Y Range"])
    XB = B.Data.index
    YB = B.Data.as_matrix(["Y Range"])
    Interpolation = interpolate.interp1d(XB, YB)
    RecastYB = Interpolation(XA)
    Subtraction = YA - RecastYB
    return FormatData(XA, Subtraction)

def InterpolateDFColumn(DataFrameA, DataFrameB):
    """ Interpolates the columns of DataFrameA to the
        index of DataFrameB.

        First match the sizes of the two dataframes,
        then interpolates the y of Frame A into x of Frame B.
    """
    CutDFA, CutDFB = NT.MatchDataFrames(DataFrameA, DataFrameB)
    NewDF = pd.DataFrame()         # New DF to store the interpolated A values
    NewX = CutDFB.index            # The values of x we will interpolate to
    OldX = CutDFA.index            # The old values of x
    for Key in CutDFA:             # Loop over columns of DataFrameA
        Interpolant = interpolate.interp1d(OldX, np.array(CutDFA[Key]))
        NewDF[Key] = Interpolant(NewX)
    NewDF.index = NewX
    return NewDF

def UpdateDictionary(OldDictionary, NewValues):
    """ Will loop over keys in new dictionary and set
    them to the old dictionary
    """
    for Key in NewValues:
        OldDictionary[Key] = NewValues[Key]
    return OldDictionary

def ConvertOrderedDict(Dictionary):
    Keys = [Key for Key in Dictionary]
    Items = [Dictionary[Key] for Key in Dictionary]
    return OrderedDict(zip(Keys, Items))

###################################################################################################

""" Fitting functions & Analysis """

def IntegrateColumn(DataFrame, Column="Y Range", Range=None):
    """ Return the integral of a column
        The range can be specified by a 2-tuple list by
        supplying the indices to begin and end.
    """
    X = np.nan_to_num(DataFrame.index)              # Ensure there are no NaNs
    Y = np.nan_to_num(np.array(DataFrame[Column]))  # in our arrays
    if Range is not None:
        X = X[Range[0]:Range[1]]
        Y = Y[Range[0]:Range[1]]
    return np.trapz(Y, X)

def FitModel(DataFrame, Model, Column="Y Range", Verbose=False):
    """ Uses an instance of the Model class to fit data contained
    in the pandas dataframe. Dataframe should have indices of the X-range
    and column "Y Range" as the Y data to be fit to

    Requires input of a DataFrame formatted using FormatData, and a reference
    to an instance of Model.

    Returns the optimised parameters, a fitting report as a dataframe,
    the fitted curves for easy plotting, and the covariance matrix.
    """
    if Verbose == True:
        print(" Curve fitting with:\t" + Model.FunctionName)
    if type(Model.BoundaryConditions) == "dict":
        Bounds = (Dict2List(Model.BoundaryConditions[0]), Dict2List(Model.BoundaryConditions[1]))
    else:
        Bounds = Model.BoundaryConditions
    if Verbose == True:
        print(" Boundary Conditions:")
        print(str(Bounds))
        print(" Initial parameters:")
    ParameterNames = [Key for Key in Model.Variables.keys()]
    OptimisedParameters, CovarianceMatrix = curve_fit(Model.Function,
                                                      np.array(DataFrame.index, dtype=float),
                                                      DataFrame[Column].values,
                                                      UnpackDict(**Model.Variables),
                                                      bounds=Bounds,
                                                      method="trf")
    ParameterReport = pd.DataFrame(data=list(zip(OptimisedParameters,
                                            CovarianceMatrix.diagonal(),
                                            np.sqrt(CovarianceMatrix.diagonal()))),
                                   columns=["Values", "Variance", "Std. Dev."],
                                   index=ParameterNames)
    ModelFit = Model.Function(np.array(DataFrame.index, dtype=float),
                              *OptimisedParameters)
    FittedCurves = pd.DataFrame(data=list(zip(DataFrame[Column], ModelFit)),
                                         columns=["Data", "Model Fit"],
                                         index=np.array(DataFrame.index, dtype=float))
    FittedCurves["Residuals"] = DataFrame[Column] - ModelFit
    if Verbose == True:
        print(" ------------------------------------------------------")
        print(" Parameter Report:")
        print(ParameterReport)
        CheckCovariance(CovarianceMatrix)
        print(" ------------------------------------------------------")
        print(" Mean-signed-error:\t" + str(np.average(np.array(FittedCurves["Residuals"]))))
        print(" RMS:\t" + str(np.average(np.square(np.array(FittedCurves["Residuals"])))))
    return OptimisedParameters, ParameterReport, FittedCurves, CovarianceMatrix

def LeastSqFit(ResidualFunction, InitialParameters, DataFrame, Column="Y Range"):
    """ Function for using scipy.optimize.leastsq to generically fit
        data, rather than use curve_fit
    """
    XData = DataFrame.index
    try:
        YData = np.array(DataFrame[Column])
    except KeyError:
        print(" No column named " + Column + " in DataFrame.")
        exit()
    return leastsq(ResidualFunction, InitialParameters, args=(XData, YData))

def ScipyLinearRegression(DataFrame, Column=None, Intercept=False, Robust=False, Margin=0.1,
                          Interface="plotly", Labels=None, Bootstrap=False):
    """ New version of a linear regression routine, coded as a wrapper
        for the least squares function from SciPy.

        Taken from the SciPy cookbook:
        http://scipy-cookbook.readthedocs.io/items/robust_regression.html

        Has the option of doing a robust regression, where outliers
        are assessed as well during the linear regression.

        FitResults is an object, containing the results of the fits.
    """
    if Intercept is False:
        def LinearResidualFunction(Parameters, XData, YData):
            """ Residual function where the objectives of the fit are
                contained in Parameters as a list or nparray

                Intercept is set to zero
            """
            return (Parameters[0] * XData) - YData
    elif Intercept is True:
        def LinearResidualFunction(Parameters, XData, YData):
            """ Residual function where the objectives of the fit are
                contained in Parameters as a list or nparray

                Intercept is fit
            """
            return (Parameters[0] * XData + Parameters[1]) - YData

    if Column == None:
        Y = np.array(DataFrame[DataFrame.keys()[0]].values)     # if none specified, do regression on the first column
    else:
        Y = np.array(DataFrame[Column].values)                  # otherwise use the specified column
    X = np.array(DataFrame.index.values)
    if Robust == False:
        LossFunction = "linear"
        Margin = None                                    # if it's normal linear regression, ignore outlier
    elif Robust == True:
        LossFunction = "soft_l1"                         # use

    InitialGuess = np.zeros(2)                            # initial guess for linear fit are ones
    FitResults = least_squares(fun=LinearResidualFunction,
                               x0=InitialGuess,
                               loss=LossFunction,
                               f_scale=Margin,
                               args=(X, Y))
    print("--------------------------------------------")
    print(" Optimal parameters:\t" + str(FitResults.x))
    FittedY = Linear(X, *FitResults.x)
    SSResiduals = np.sum(np.square(FittedY - Y))
    SSY = np.sum(np.square(Y - np.average(Y)))
    print(" R**2 value:\t" + str(1. - SSResiduals / SSY))
    print("--------------------------------------------")
    if Bootstrap == True:
        Parameters = []
        NTrials = 100
        for Trial in xrange(NTrials):
            SimulatedY = AddNoise(FittedY)
            BootstrapResults = least_squares(fun=LinearResidualFunction,
                                             x0=InitialGuess,
                                             loss=LossFunction,
                                             f_scale=Margin,
                                             args=(X, SimulatedY))
            Parameters.append(BootstrapResults.x)
        print(" Bootstrap analysis results:")
        print(" Averages:\t" + str(np.average(Parameters, axis=0)))
        print(" Sigma:\t" + str(np.std(Parameters, axis=0)))
        FitResults.BootstrapAverage = np.average(Parameters, axis=0)
        FitResults.BootstrapStdDev = np.std(Parameters, axis=0)
    print("--------------------------------------------")
    DataFrame = pd.DataFrame(data=list(zip(Y, FittedY)),
                             index=X,
                             columns=["Y Range", "Linear Regression"])
    PlotData(DataFrame,
             Interface=Interface,
             Labels=Labels
             )
    return DataFrame, FitResults

def LinearRegression(DataFrame, Columns=None, Labels=None):
    RegressionReport = dict()
    LinearModel = Model("Linear regression model")       # Set up regression model
    if Columns is None:
        Columns = list(DataFrame.keys())            # Default is regress all columns
    with NT.suppress_stdout():
        # Supress output
        LinearModel.SetFunction(Linear)
        for Data in Columns:
            if "-Regression" in Data:         # No need to regress more than once!
                pass
            else:
                Opt, Report, Fits, Cov = FitModel(DataFrame, Model=LinearModel, Column=Data, Verbose=False)
                DataFrame[Data + "-Regression"] = Linear(DataFrame.index, *Opt)    # Plot the curves
                RegressionReport[Data] = Report
                if (Data + "-Regression") not in Columns:
                    Columns.append(Data + "-Regression")
    fig, ax = PlotData(DataFrame=DataFrame, Columns=Columns, Labels=Labels)
    return RegressionReport

def CheckCovariance(CovarianceMatrix, Threshold=1E-2):
    OffDiagonalElements = NT.CheckOffDiagonal(CovarianceMatrix)
    ThresholdedElements = [Element for Element in OffDiagonalElements if Element >= Threshold]
    if len(ThresholdedElements) > 0:
        print(str(len(ThresholdedElements)) + " elements of covariance matrix larger than\t" + str(Threshold))
        print(" Check your privilege!")

def ProgressionFitting(PeakNumbers, PeakReport, Offset=0):
    """ Function for fitting a second order polynomial to a set of peaks
        in some hope of finding out what the shit Ethane is...

        Takes PeakNumbers as a list or array of integers corresponding to the
        peak numbers we want to assign to a progression, and PeakReport which
        holds the locations, indices and intensities for the peaks.

        Offset will shift the M counter (arbitrary J effectively).
    """
    Locations = []
    NPeaks = len(PeakNumbers)
    for Peak in PeakNumbers:
        Locations.append(PeakReport.iloc[Peak]["X Value"])              # Find the actual peak values

    """ Set up the model """
    PolynomialModel = Model("2nd order polynomial fit")
    PolynomialModel.SetFunction(SecondOrderPolynomial)                  # set fitting function as quadratic

    ProgressionObject = Spectrum(Reference="Progression")
    ProgressionObject.Data = pd.DataFrame(data=Locations, index=np.arange(Offset, Offset+NPeaks), columns=["Y Range"])
    ProgressionObject.Labels = {"X Label": "M","Y Label": "Energy","Title": "Progression Fit"}
    ProgressionObject.Fit(PolynomialModel, Verbose=True,)
    return ProgressionObject

def PeakFinding(DataFrame, Column=None, Threshold=0.3, MinimumDistance=30.):
    """ Routine that will sniff out peaks in a spectrum, and fit them with Gaussian functions
    I have no idea how well this will work, but it feels like it's gonna get pretty
    complicated pretty quickly
    """
    if Column is None:
        Column = DataFrame.keys()[0]         # if nothing supplied, use the first column
        print("No column specified, using " + Column)
    PeakIndices = peakutils.indexes(DataFrame[Column], thres=Threshold, min_dist=MinimumDistance)
    NPeaks = len(PeakIndices)
    print(" Found \t" + str(NPeaks) + "\t peaks.")
    StickSpectrum = np.zeros((len(DataFrame[Column])), dtype=float)
    Intensities = []
    for Index in PeakIndices:
        StickSpectrum[Index] = 1.
    DataFrame["Stick Spectrum"] = StickSpectrum
    return PeakIndices

def VMICalibration(Spectrum):
    """ Procedure for VMI calibration. This way we don't have
        to redo it over and over again.

        1. Detect peaks in loaded spectrum
        2. Assign peak indices to atom levels
    """
    Threshold = float(raw_input(" Threshold for peak detection:"))
    Spectrum.DetectPeaks(Threshold=Threshold)
    Spectrum.PlotLabels()
    #from bokeh.io import output_notebook
    #output_notebook()
    Spectrum.PlotAll(Interface="plotly",
                     Columns=["Y Range"],
                     Labels={"X Label": "Pixel speed",
                             "Y Label": "P(s)",
                             "Title": "Uncalibrated speed distribution"
                             },
                     PlotTypes={"Y Range": "line"},
                     )
    NLines = int(raw_input(" Number of peaks used for calibration?"))
    PeakAssignments = ["3P", "5P", "3S", "5S"]
    PeakDict = OrderedDict()
    for n in xrange(NLines):
        PeakDict[PeakAssignments[n]] = int(raw_input(" Index for peak:\t" + PeakAssignments[n]))
    CalibrationData = FormatData(X=sorted(Dict2List(PeakDict)), Y=sorted(Dict2List(OxygenAtomSpeed)))
    LinearRegression = Model("Linear Regression")
    LinearRegression.SetFunction(Linear)
    LinearRegression.SetVariables({"Gradient": 16.,
                                   "Offset": 10.})
    try:
        popt, report, fits, pcov = FitModel(CalibrationData, LinearRegression)
        Labels = {"Title": "Calibration for:    " + Spectrum.Reference,
              "X Label": "Pixel speed",
              "Y Label": "Oxygen atom speed",
              }
        PlotData(fits, Labels)
    except TypeError:
        print(" Could not automatically fit data. Please try it manually. ")
        report = 0.
    finally:
        return report, fits

def DataFrameDistances(DataFrame, PeaksRange=[0,5]):
    DistancesDataFrame = pd.DataFrame()
    for PeakNumber in range(PeaksRange[0], PeaksRange[1]):       # Loop over all indices
        PeakPosition = DataFrame.index[PeakNumber]
        DistancesDataFrame[PeakNumber] = np.abs(DataFrame.index - PeakPosition)
    return DistancesDataFrame

def AddNoise(Y, DampFactor=0.5):
    """
    Function to add random noise to data; this will be used for Bootstrap
    analysis. Superior to the previously written version because
    it won't add noise to where it shouldn't be.

    Should work for any arbitrary size of data

    A damp factor is included to soften the noise if it's too noisy.
    """
    NewData = np.zeros((len(Y)), dtype=float)
    Iterator = np.nditer(Y, flags=['f_index'])      # loop over elements
    while not Iterator.finished:
        Delta = np.random.rand() >= 0.5             # flip a coin to determine
        if Delta == True:                           # the sign of the noise
            Multiplier = 1. * DampFactor
        else:
            Multiplier = -1. * DampFactor
        NewData[Iterator.index] = Y[Iterator.index] + (Y[Iterator.index] * np.random.rand() * Multiplier)
        Iterator.iternext()                         # Next iteration of loop
    return NewData

def HistogramBin(Array, Bins="auto", KDE=True, Plot=False):
    Histogram, Bins = np.histogram(Array,
                                   bins=Bins,
                                   density=KDE,
                                   )
    DataFrame = pd.DataFrame(data=Histogram, index=Bins[:-1], columns=["Histogram"])
    if Plot is True:
        PlotData(DataFrame, PlotTypes={"Histogram": "bar"})
    else:
        pass
    return DataFrame

def AnisotropyAveraging(Spectrum, PixelRange=[165,195]):
    """ Calculates the average anisotropy value for a range of pixels
        in an angular distribution.
    """
    BetaValues = []                                   # list that stores the beta values
    for Pixel in range(PixelRange[0], PixelRange[1]):                   # loop over pixels
        PixelValue = Spectrum.Data.index[NT.find_nearest(Spectrum.Data.index, Pixel)]
        BetaValues.append(Spectrum.Data["Beta"][PixelValue])       # find the beta value at a given pixel
    BetaValues = np.array(BetaValues)
    return np.average(BetaValues, axis=0), np.std(BetaValues, axis=0)

###################################################################################################

""" Commonly used base functions """

def BaseGaussian(x, x0):
    """ The most vanilla of Gaussians, wrote it when I was debugging
    """
    return np.exp(-np.square(x-x0))

def ScottGaussian(x, Amplitude, Centre, Width):
    """ This version of the Gaussian I used for the convolution fitting
        Notably, the width is given as the full width, not half!
    """
    return Amplitude * np.exp(-(Centre - x)**2. / Width**2.)

def GaussianFunction(x, Amplitude, Centre, Width):
    return Amplitude * (1 / (Width * np.sqrt(2 * np.pi))) * np.exp(-np.square(x - Centre) / (2 * Width**2))

def BoltzmannFunction(x, Amplitude, Temperature):
    #return Amplitude * np.exp(1) / (kcm * Temperature) * x * np.exp(- x / (kcm * Temperature))
    return Amplitude * np.sqrt(1 / (2 * np.pi * kcm * Temperature)**3) * 4 * np.pi * x * np.exp(-(x) / (kcm * Temperature))

def SecondOrderPolynomial(x, A, B, C):
    return (A * x**2.) + (B * x) + C

def Linear(x, Gradient, Offset):
    return x * Gradient + Offset

def PriorFunction(x, A, i, j):
    return A * (x**i) * (1. - x)**j

def SingleExponential(x, Amplitude, Gradient, Offset):
    return Amplitude * np.exp(Gradient * (x + Offset))

@jit
def ConvolveArrays(A, B, method="new"):
    """ Function I wrote to compute the convolution of two arrays.
    The new method is written as a manual convolution calculated by
    taking the inverse Fourier transform of the Fourier product of the
    two arrays.

    Returns only the real values of the transform.

    Old function uses the scipy function, which I had issues with the
    length of array returned.

    Requires input of two same-dimension arrays A and B, and an optional
    1-D array that holds the X dimensions.
    """
    if method == "old":
        BinSize = len(A)
        ConvolutionResult = signal.fftconvolve(A, B, mode="full")
        ReshapedConvolution = ConvolutionResult[0:BinSize]
        return ReshapedConvolution                                 # Return same length as input X
    elif method == "new":
        FTA = fftpack.fft(A)
        FTB = fftpack.fft(B)
        FTProduct = FTA * FTB
        return np.real(fftpack.ifft(FTProduct))

###################################################################################################

""" Commonly used model objects """

###################################################################################################

""" File I/O and plotting functions """

def LoadSpectrum(File, CalConstant=1.):
    """ Function for generating a pandas dataframe from a csv,
    involves smart sniffing of delimiter, and returns a dataframe where
    the index is the X Range of the spectrum.

    Assumes two-column data.
    """
    Delimiter = NT.DetectDelimiter(File)                          # detect what delimiter
    if NT.DetectHeader(File) != True:                             # see if header is there
        DataFrame = pd.read_csv(File, delimiter=Delimiter,
                                header=None, names=["Y Range"],
                                index_col=0)                      # read file from csv
    else:                                                         # if header is present
        DataFrame = pd.read_csv(File, delimiter=Delimiter,
                                index_col=0)
    DataFrame.set_index(DataFrame.index * CalConstant, inplace=True) # Modify index
    DataFrame = DataFrame.dropna(axis=0)                          # removes all NaN values
    return DataFrame

def DatabaseSpectrum(Database, Reference):
    Data = NT.LoadReference(Database, Reference)
    Filename = raw_input("Please specify which spectrum to load.")
    SelectedData = Data[Filename]
    DataFrame = pd.DataFrame(data=SelectedData[:,2],
                             index=SelectedData[:,1],
                             columns=["Y Range"])
    return DataFrame

def TDLWavenumber(DataFrame):
    Wavenumber = 2e7 / DataFrame.index + 1e7 / (1064.464)
    DataFrame.index = Wavenumber

def GenerateJComb(DataFrame, TransitionEnergies, Offset=1.3, Teeth=0.2, SelectJ=10):
    """ Function for generating a rotational comb spectrum as an annotation
    Not elegant, but it works!
    The required input:
    DataFrame - pandas dataframe containing X axis of target spectrum
    TransitionEnergies - A list of transition energies in some (e.g. J) order
    Offset - The height of the comb
    Teeth - The length of the comb's teeth (lol)
    SelectJ - Integer for selecting out multiples of J so it's not insane

    Sets the input dataframe key "Comb" with the comb spectrum
    """
    # Find the indices closest to where we can put our comb teeth on
    #Indices = [NT.find_nearest(DataFrame.index, Energy) for Energy in TransitionEnergies if TransitionEnergies[Energy] % SelectJ == 0]
    Indices = []
    for index, Energy in enumerate(TransitionEnergies):
        if index % SelectJ == 0:
            Indices.append(NT.find_nearest(DataFrame.index, Energy))
    Comb = [Offset for value in DataFrame.index]
    print(" Adding\t" + str(len(Indices)) + "\t teeth to the comb.")
    for index in Indices:
        Comb[index] = Comb[index] - Teeth
    DataFrame["Comb"] = Comb

def PickleSpectra(DataBase):
	PickleDictionary = dict()
	for Instance in Spectrum.instances:
		try:
			PickleDictionary[Instance.Reference] = Instance
		except AttributeError:
			print(" No reference for " + Instance + ", skipping.")
			pass
	NT.SaveObject(PickleDictionary, DataBase)

def SetupPyplotTypes(Key, Scaling=None, Colours=None, PlotType=None):
    """ Programatically select the kind of plots we can
        produce. Use PlotType as input for the kind of
        plot, and this function will return the appropriate pyplot
        function as well as a dictionary of settings

        Colours is a dictionary containing the colour of each bar and
        scatter plot.

        Key is the name of the column in the DataFrame.

        Scaling is a tuple which dynamically scales the size of markers
        according to the number of data points.
    """
    Wrappers = {"scatter": plt.scatter,         # dictionary of pyplot functions
                "bar": plt.bar,
                "line": plt.plot}
    DefaultSettings = {"scatter": {"s": 150. * Scaling,  # settings for each
                                   "alpha": 0.7,        # kind of plot
                                   "label": Key,
                                   "c": Colours[Key]
                                  },
                        "line": {"alpha": 0.8,
                                 "antialiased": True,
                                 "label": Key},
                        "bar": {"alpha": 0.8,
                                "width": 0.5,
                                "align": "center",
                                "color": Colours[Key],
                                "label": Key}
                       }
    """ Set up the initial default settings before modifying """
    if PlotType == None:                        # default plot type is scatter
        Settings = DefaultSettings["scatter"]
    else:
        Settings = DefaultSettings[PlotType]
    try:
        print(Settings[PlotType]["c"])
    except KeyError:
        pass
    return Wrappers[PlotType], Settings

def PlotData(DataFrame, Labels=None, Columns=None, Legend=True, Interface="pyplot", PlotTypes=None, Annotations=None):
    """ A themed data plotting routine. Will use either matplotlib or
    bokeh to plot everything in an input dataframe, where the index is
    the X axis.

    PlotTypes is a dictionary input with the same keys as the DataFrame
    columns.

    Labels is a dictionary containing the axis plot settings.
    """
    if Columns is None:
        Columns = DataFrame.keys()
    Headers = DataFrame.keys()                  # Get the column heads
    if Interface == "pyplot":                                 # Use matplotlib library
        fig = plt.figure(figsize=(12,8))
        plt.margins(0.2)                                # Margins from the edges of the plot
        #plt.axes(frameon=True)
        ax = fig.gca()

        """ Set up the kinds of plots that will be displayed.
            Plots is a dictionary that holds the kind of plot for each
            column in DataFrame.
        """
        Plots = dict()
        for Key in Columns:
            if (type(Key) is str) == True:       # check if header is string, skip if not
                Exists = NT.CheckString(Key, ["Model", "Regression", "Fit", "Smoothed"])
                if Exists is True:
                    Plots[Key] = "line"        # if the data is actually a fitted model
                elif Exists is False:
                    Plots[Key] = "scatter"     # Initialise default plots to scatter

        """ If there are specific plot types requested,
            update the dictionary.
        """
        if PlotTypes is not None:           # Only update if there are
            for Key in PlotTypes:          # specified plot types
                Plots[Key] = PlotTypes[Key]
        else:
            pass

        """ Set up the plot labels for the axes. """
        if Labels is not None:
            try:                                              # Unpack whatever we can from Labels
                plt.xlabel(Labels["X Label"], fontsize=18.)
                plt.ylabel(Labels["Y Label"], fontsize=18.)
                plt.title(Labels["Title"])
                plt.xlim(Labels["X Limits"])
                plt.ylim(Labels["Y Limits"])
            except KeyError:                    # Will ignore whatever isn't given in dictionary
                pass

        """ Set up the colours used for scatter and bar plots """
        Colours = dict()
        ColourCounts = 0
        for Key in Plots:                                # Count the number of bar/scatter plots
            Colours[Key] = "blue"
            if Plots[Key] is "scatter" or "bar":
                ColourCounts = ColourCounts + 1
        ColourMap = cm.Accent(np.linspace(0, 1, ColourCounts))   # Interpolate colours
        for Index, Key in enumerate(Plots):                      # Assign colour to plot
            Colours[Key] = ColourMap[Index]

        Scaling = 1. / ((DataFrame.shape[0]))**(1. / 12.)   # Scale the marker sizes to number of elements
        """ Loop over each column in the DataFrame. """
        for Key in Columns:
            """ Determine what kind of plot is needed, and the
                settings to go with it.
            """
            PlotFunction, PlotSettings = SetupPyplotTypes(Key=Key,
                                                          Scaling=Scaling,
                                                          PlotType=Plots[Key],
                                                          Colours=Colours,
                                                          )
            PlotFunction(DataFrame.index, DataFrame[Key], **PlotSettings)
        if Legend is True:
            plt.legend(ncol=1, loc=0)
        ax.grid()
        #ax.set_autoscale(enable=True, tight=True)
        plt.show()
    elif Interface == "bokeh":                                # Use bokeh library
        NCols = len(DataFrame.columns)                            # Get the number of columns
        if NCols <= 2:
            Colours = ["blue", "red"]
        else:
            Colours = brewer["Spectral"][NCols]                   # Set colours depending on how many
        tools = "pan, wheel_zoom, box_zoom, reset, resize, hover"
        ax = None                                             # Not used with Bokeh
        if Labels != None:
            try:                                              # Unpack whatever we can from Labels
                XLabel = Labels["X Label"]
                YLabel = Labels["Y Label"]
                Title = Labels["Title"]
                XRange = Labels["X Limits"]
                YRange = Labels["Y Limits"]
                fig = figure(width=700, height=400,                        # set up the labels
                              x_axis_label=XLabel, y_axis_label=YLabel,
                              title=Title, tools=tools,
                              x_range=XRange, y_range=YRange)
            except KeyError:
                print(" Not using labels")
                pass
        else:
            fig = figure(width=700, height=400, tools=tools)               # if we have no labels
        fig.background_fill_color="beige"
        for index, Data in enumerate(DataFrame):
            fig.scatter(x=DataFrame.index, y=DataFrame[Data],
                      line_width=2, color=Colours[index],
                      legend=Headers[index])
            fig.line(x=DataFrame.index, y=DataFrame[Data],
                      line_width=2, color=Colours[index],
                      legend=Headers[index])
        show(fig)
    elif Interface == "plotly":                                # Plotly interface
        import PlottingTools as PT
        PT.XYPlot(DataFrame=DataFrame,
                  Columns=Columns,
                  CustomPlotTypes=PlotTypes,
                  Labels=Labels,
                  Annotations=Annotations
                  )
        fig, ax = (None, None)
    return fig, ax

def ToggleD3():
    import mpld3
    if D3State is False:
        mpld3.enable_notebook()
    elif D3State is True:
        mpld3.disable_notebook()
