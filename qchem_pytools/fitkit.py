
from scipy.optimize import leastsq, curve_fit
from qchem_pytools import figure_settings
import numpy as np
import pandas as pd

class leastsq_object:
    def __init__(self):
        print("Usage: specify model function as leastsq_object.model")
        print("The model should have parameters as an n-tuple, followed by x.")
        self.model = None

    def setup(self, model=None, initial=None):
        if model is not None:
            """ define some stock functions """
            if model == "linear":
                self.model = lambda x, param: param[0] * x + param[1]
            elif model == "quadratic":
                self.model = lambda x, param: param[0] * x**2. + param[1] * x + param[2]
            elif model == "gaussian":
                self.model = lambda x, param: param[0] * np.exp(-(x - param[1])**2. / (2. * param[2]**2.))
            elif modell == "cbs":
                self.model = lambda x, param: param[0] + param[1] * x**(-3.)
            else:
                self.model = model
            """ initialise error function """
            self.error_function = lambda param, x, y: self.model(x, param) - y
            if initial is not None:
                self.initial = initial
            else:
                self.initial = np.zeros(6)
        elif model is None:
            print("Model function undefined, provide as argument.")

    def optimise(self, x ,y):
        if self.model is not None:
            optparam, success = leastsq(
                func=self.error_function,
                x0=self.initial,
                args=(x, y)
            )
            if success >= 1 and success <= 4:
                print("Optimisation succeeded.")
                self.modeled = self.model(x, optparam)
                self.residuals = self.modeled - y
                self.sos = np.sum(self.residuals**2.)
                self.optparam = optparam
                self.success = True
                print("Final parameters:\t" + str(optparam))
                print("Sum of squares:\t" + str(self.sos))
                self.dataframe = pd.DataFrame(
                    data=[self.modeled, y, self.residuals],
                    index=x,
                    columns=["Modelled", "Actual", "Residuals"]
                    )
            else:
                print("Optimisation failed.")
        else:
            print("Model function undefined, call setup method")

    def plot(self):
        subplot_grid = figure_settings.GenerateSubPlotObject(
                ColumnNames=["Fit", "Residuals"],
                NRows=2,
                NCols=1
            )
        actual_trace = figure_settings.DefaultScatterObject(
                X=self.dataframe.index,
                Y=self.dataframe["Actual"],
                Name="Actual",
            )
        modelled_trace = figure_settings.DefaultScatterObject(
                X=self.dataframe.index,
                Y=self.dataframe["Modelled"],
                Name="Modelled",
                Connected=True
            )

        residual_trace = figure_settings.DefaultScatterObject(
                X=self.dataframe.index,
                Y=self.dataframe["Residuals"],
                Name="Residuals"
            )
        subplot_grid.append_trace([actual_trace, modelled_trace], 1, 1)
        subplot_grid.append_trace(residual_trace, 2, 1)
        figure_settings.iplot(subplot_grid)

class curvefit_object:
    def __init__(self):
        print("Usage: specify model function as curvefit_object.model")
        print("The model should have parameters as an n-tuple, followed by x.")
        self.model = None
        self.bounds = {"upper": [], "lower": []}

    def setup(self, model=None, initial=None, bounds=None):
        if model is not None:
            """ define some stock functions here """
            if model == "linear":
                self.model = lambda x, param: param[0] * x + param[1]
            elif model == "quadratic":
                self.model = lambda x, param: param[0] * x**2. + param[1] * x + param[2]
            elif model == "gaussian":
                self.model = lambda x, param: param[0] * np.exp(-(x - param[1])**2. / (2. * param[2]**2.))
            elif modell == "cbs":
                self.model = lambda x, param: param[0] + param[1] * x**(-3.)
            else:
                self.model = model
        else:
            raise("ModelError! No model specified as argument!")

    def fit(self, x, y):
        if self.model is not None:
            optparam, covar = curve_fit(
                f=self.model,
                xdata=x,
                ydata=y,
                p0=self.initial
            )
            colours = figure_settings.GenerateColours(NumPlots=2)
            actual_trace = figure_settings.DefaultScatterObject(
                X=x,
                Y=y,
                Name="Data",
                Colour=colours[0]
            )
            model_trace = figure_settings.DefaultScatterObject(
                X=x,
                Y=self.model(x, *optparam),
                Name="Fit",
                Colour=colours[1]
            )
            layout = figure_settings.DefaultLayoutSettings()
            figure_settings.iplot(
                {
                    "data": [actual_trace, model_trace],
                    "layout": layout
                }
            )
        else:
            raise("ModelError! No model defined!")
