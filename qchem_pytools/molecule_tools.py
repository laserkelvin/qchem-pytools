
from scipy import constants
import numpy as np
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from .atom_data import masses, radii

# Since I can't get PyBel to work properly, gonna try and
# write a few of my own routines to get the ball rolling

############################### Classes ###############################

# Molecule class, a dictionary of atoms. Might be a bit clunky to 
# reference...
class Molecule:
    # Initialises by adding instances of atoms class to the dictionary
    # from XYZ format
    def __init__(self, InputXYZ=None):
        """ Initialises an instance of Molecule
        Keyword arguments:
        InputXYZ - The XYZ coordinates as read in using the function ReadXYZ
        Attributes:
        AtomScalingFactor - Used for scaling the size of atoms in Show() method
        COMCoordinates - Boolean flag for the coordinates are in COM frame
        CalculatedCOM - Boolean flag for whether or not COM has been calculated 
        Atoms - Dictionary that holds instances of Atoms class 
        COM - np.ndarray holding COM coordinates
        
        Methods:
        CalculateCOM() - Returns the center-of-mass coordinates, and flags
                         the COM coordinates to be calculated.
        Shift2COM() - Shifts the coordinates of all the atoms from arbitrary
                      laboratory frame to the centre-of-mass frame.
        ExportXYZ() - Returns the XYZ coordinates in a way that's copy-pastable
                      to wxMacMolPlt or some electronic structure program
        Show() - Plots the atoms of the molecule in XYZ coordinates using the
                 matplotlib scatter library. Additional things that are plotted:
                 If the COM has been calculated:
                 it will also show this as a purple sphere. 
                 If the PrincipalAxes have been calculated:
                 it will show the colour coded axes - RGB for elements 0,1,2
        GenerateBonds() - Returns connectivity between atoms based on a threshold
                          minimum distance
        ChemView() - Uses the ChemViewer module in IPython Notebook to display a
                     3D model of the molecule. Requires GenerateBonds() to pipe 
                     connectivity.
        CalculateInertiaMatrix() - Returns the moments of inertia (diagonal) and
                                   products of inertia (off-diagonal) matrix in
                                   SI units (kg m**2). This method assumes the
                                   molecule is already in a COM frame by Shift2COM.
        PrincipalMoments() - Returns the PrincipalMoments and PrincipalAxes by
                             diagonalising the inertia matrix. PrincipalMoments are
                             reported in 1/cm.
        """
        self.AtomScalingFactor = 600.                           # Used for method Show()
        self.COMCoordinates = False                             # Initial state is not in COM coordinates
        self.CalculatedCOM = False                              # Whether or not COM has been calculated already
        self.Atoms = {}                                         # Dictionary holding instances of Atom class
        self.COM = np.array([0., 0., 0.])                       # COM of molecule
        if InputXYZ is not None:
            """ If no input is specified, don't read in """
            self.Atoms = MoleculeFromXYZ(InputXYZ)
            self.NAtoms = len(self.Atoms)
    # Function to calculate the centre of mass for a given molecule in XYZ coordinates
    def CalculateCOM(self):
        CumSum = 0.                                             # Cumulative sum of m_i * sum(r_i)
        TotalMass = 0.                                          # Total mass of molecule
        for Atom in self.Atoms:
            CumSum = CumSum + self.Atoms[Atom].GetMass() * self.Atoms[Atom].GetCoordinates()
            TotalMass = TotalMass + self.Atoms[Atom].GetMass()
        #for AtomNumber in range(self.NAtoms):               # Loop over all atoms
        #    CumSum = CumSum + self.Atoms[str(AtomNumber)].Mass * self.Atoms[str(AtomNumber)].Coordinates()
        #    TotalMass = TotalMass + self.Atoms[str(AtomNumber)].Mass
        self.CalculatedCOM = True                               # Flag COM as having been calculated for plotting
        self.COM = (1. / TotalMass) * CumSum
        return self.COM
    # Function to shift the coordinates to centre of mass frame (I think it works)
    # by making the COM zero
    def Shift2COM(self):
        if self.COMCoordinates == True:
            print(" Already in COM frame")
            pass
        else:
            self.CalculateCOM()
            for Atom in self.Atoms:
                self.Atoms[Atom].Coordinates = self.Atoms[Atom].GetCoordinates() + (-self.COM)
            self.COMCoordinates = True
    # Function that will format xyz
    def ExportXYZ(self):
        for Atom in self.Atoms:
            print(self.Atoms[Atom].GetSymbol() + "\t" + str(self.Atoms[Atom].GetCoordinates()))

    # Function to plot up the molecule using an xyz matplotlib plot
    def Show(self):
        HydrogenRadius = 53.                                    # in picometres
        #AtomicRadii = {"H": 53. / HydrogenRadius,
        #               "C": 67. / HydrogenRadius,
        #               "O": 48. / HydrogenRadius,
        #               "N": 56. / HydrogenRadius,
        #               "Na": 190. / HydrogenRadius,
        #               "Mg": 145. / HydrogenRadius,
        #               "Al": 118. / HydrogenRadius,
        #               "Si": 111. / HydrogenRadius,
        #               "P": 98. / HydrogenRadius,
        #               "S": 88. / HydrogenRadius,
        #               "K": 243. / HydrogenRadius,
        #               "Ca": 194. / HydrogenRadius,
        #               "Sc": 184. / HydrogenRadius,
        #               "Ti": 176. / HydrogenRadius,
        #               "V": 171. / HydrogenRadius,
        #               "Cr": 166. / HydrogenRadius,
        #               "Mn": 161. / HydrogenRadius,
        #               "Fe": 156. / HydrogenRadius,
        #               "Co": 152. / HydrogenRadius,
        #               "Ni": 149. / HydrogenRadius,
        #               "Cu": 145. / HydrogenRadius,
        #               "Zn": 142. / HydrogenRadius,
        #               "COM": 50. / HydrogenRadius,
        #               "X": 50. / HydrogenRadius}
        AtomicColours = {"H": "#FFFFFF",                          # CPK colours
                         "C": "#909090",
                         "N": "#3050F8",
                         "O": "#FF0D0D",
                         "Na": "#AB5CF2",
                         "Mg": "#8AFF00",
                         "Al": "#BFA6A6",
                         "Si": "#F0C8A0",
                         "P": "#FF8000",
                         "K": "#8F40D4",
                         "Ca": "#3DFF00",
                         "Sc": "#E6E6E6",
                         "Ti": "#BFC2C7",
                         "V": "#A6A6AB",
                         "Cr": "#8A99C7",
                         "Mn": "#9C7AC7",
                         "Fe": "#E06633",
                         "Co": "#F090A0",
                         "Ni": "#50D050",
                         "Cu": "#C88033",
                         "Zn": "#7D80B0",
                         "COM": "green",
                         "X": "purple"}
        NAtoms = len(self.Atoms)
        Colors = []
        if self.CalculatedCOM == False:
            X = np.zeros((NAtoms), dtype=float)                 # Arrays for holding xyz coordinates
            Y = np.zeros((NAtoms), dtype=float)
            Z = np.zeros((NAtoms), dtype=float)
            Size = np.zeros((NAtoms), dtype=float)              # Array for holding size of atoms
        else:
            X = np.zeros((NAtoms + 1), dtype=float)             # Arrays for holding xyz coordinates
            Y = np.zeros((NAtoms + 1), dtype=float)             # one more element for COM point
            Z = np.zeros((NAtoms + 1), dtype=float)
            Size = np.zeros((NAtoms + 1), dtype=float)          # Size of the atoms
        for Index, AtomNumber in enumerate(self.Atoms):                        # Loop over all atoms
            X[Index], Y[Index], Z[Index] = self.Atoms[AtomNumber].GetCoordinates()
            AtomicSymbol = self.Atoms[AtomNumber].GetSymbol()
            Colors.append(AtomicColours[AtomicSymbol])          # work out the colour for atom
            Size[Index] = (float(radii[AtomicSymbol]) / HydrogenRadius) * self.AtomScalingFactor
        if self.CalculatedCOM == True:                          # If we calculated COM before plot it too
            X[NAtoms], Y[NAtoms], Z[NAtoms] = self.COM
            Size[NAtoms] = AtomicRadii["COM"] * self.AtomScalingFactor
            Colors.append(AtomicColours["COM"])
        fig = plt.figure()                                      # Use matplotlib to plot atoms in 3D scatter
        ax = plt.axes(projection = "3d")
        ax.w_xaxis.gridlines.set_lw(5.0)                        # This makes the grid lines thicker apparently
        ax.w_yaxis.gridlines.set_lw(5.0)
        ax.w_zaxis.gridlines.set_lw(5.0)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.scatter(X, Y, Z, s=Size, c=Colors, depthshade=False)
        try:                                                    # if the principal axes have been calculated,
            type(self.PrincipalAxes)                            # show them on the plot too
            PX = np.zeros((3), dtype=float)                     # arrays for holding the principal axes
            PY = np.zeros((3), dtype=float)                     # components
            PZ = np.zeros((3), dtype=float)
            AxesColors = ["red", "green", "blue"]
            for AxisNumber in enumerate(self.PrincipalAxes):
                Vector = AxisNumber[1]
                ax.plot([0, Vector[0]], [0, Vector[1]], [0, Vector[2]],
                antialiased=True, color=AxesColors[AxisNumber[0]], linestyle="dashed")
        except AttributeError:
            pass

        plt.show()
    def GenerateBonds(self, Threshold=1.6):                     # Generate a list of bonds based on distance
        Bonds = []
        for Atom1 in self.Atoms:
            for Atom2 in self.Atoms:
                if Atom1 == Atom2:                  # If it's the same atom don't worry about it
                    pass
                else:
                    Distance = CalculateDistance(self.Atoms[Atom1], self.Atoms[Atom2])
                    if Distance <= Threshold:                   # If distance is less than threshold, it's a bond!
                        Bonds.append((Atom1, Atom2))
        self.Bonds = Bonds
        return Bonds
    # Routine effectively copied from down below, but adapted for use as molecule method
    # Calculates the moment of inertia matrix without assuming that it's diagonal...
    def CalculateInertiaMatrix(self):
        """
        Calculates the Inertia Tensor assuming a center-of-mass frame.
        """
        if self.COMCoordinates == False:                        # come up with warning about not being in COM frame
            print(" Warning: not in COM frame! Inertia matrix will be wrong!")
        else:
            pass                                                # no need to worry if already in COM frame
        InertiaMatrix = np.zeros((3,3), dtype=float)
        self.InertiaMatrix = InertiaMatrix                      # Zero the matrix beforehand
        for Atom in self.Atoms:                   # Loop over atoms in molecule
            Coordinates = self.Atoms[Atom].GetCoordinates() * 1e-10      # Retrieve coordinates, convert to metres
            Mass = self.Atoms[Atom].GetMass()                            # Retrieve mass of atom
            # Work out diagonal elements of the matrix
            InertiaMatrix[0,0] = InertiaMatrix[0,0] + IXX(Coordinates[1], Coordinates[2], Mass)
            InertiaMatrix[1,1] = InertiaMatrix[1,1] + IYY(Coordinates[0], Coordinates[2], Mass)
            InertiaMatrix[2,2] = InertiaMatrix[2,2] + IZZ(Coordinates[0], Coordinates[1], Mass)
            # Calculate off-diagonal elements of matrix
            InertiaMatrix[0,1] = InertiaMatrix[0,1] + IXY(Coordinates[0], Coordinates[1], Mass)
            InertiaMatrix[0,2] = InertiaMatrix[0,2] + IXZ(Coordinates[0], Coordinates[2], Mass)
            InertiaMatrix[1,2] = InertiaMatrix[1,2] + IYZ(Coordinates[1], Coordinates[2], Mass)
        # Apply sign change
        InertiaMatrix[0,1] = -InertiaMatrix[0,1]
        InertiaMatrix[0,2] = -InertiaMatrix[0,2]
        InertiaMatrix[1,2] = -InertiaMatrix[1,2]
        # Symmetrize 
        self.InertiaMatrix = (InertiaMatrix + InertiaMatrix.T) / 2.
        return self.InertiaMatrix
    # Diagonalise the moments of inerta matrix and work out the principal moments of inertia
    def PrincipalMoments(self):
        if type(self.InertiaMatrix) == None:                    # check if MOI matrix has been calculated
            self.CalculateInertiaMatrix()
        else:
            Diagonal = np.linalg.eig(self.InertiaMatrix)     # Ignore the eigenvectors
            Eigenvalues = Diagonal[0]
            Eigenvectors = Diagonal[1]
            self.PMI = PMI2ABC(Eigenvalues)                  # Return the rotational constants in 1/cm
            self.PrincipalAxes = Eigenvectors                # Return the eigenvectors for plotting
            return self.PMI
    def SumMass(self):                                          
        Mass = 0.
        for AtomNumber in range(self.NAtoms):
            Mass = Mass + self.Atoms[str(AtomNumber)].Mass
        self.Mass = Mass
        return self.Mass

# Atom class, has attributes of the xyz coordinates as well as its symbol and mass
class __Atom__:
    def __init__(self, Coordinate, Symbol):
        """ Label is atom symbol string,
            Coordinates is 3-tuple list
            
            Charges are taken from this list:
            https://en.wikipedia.org/wiki/Effective_nuclear_charge#Values
        """
        ChargeTable = {"C": 12.026,
                       "B": 9.677,
                       "H": 1.,
                       "O": 16.603,
                       "N": 14.346,
                       "He": 1.688,
                       "X": 0.0,               # Dummy atom
                      }
        self.Coordinates = np.array(Coordinate)
        self.Symbol = Symbol
        self.Mass = Symbol2Mass(self.Symbol)
        if self.Symbol in ChargeTable.keys():
            self.Charge = ChargeTable[self.Symbol]
        else:
            self.Charge = 0.
    # Returns the coordinates of an atom
    def GetCoordinates(self):
        return self.Coordinates
    
    def GetSymbol(self):
        return self.Symbol
    
    def GetCharge(self):
        return self.Charge

    def GetMass(self):
        return self.Mass

############################### I/O ###############################

# Function to read the output xyz coordinates. The below functions
# will work better if this is rotated to the inertial frame of reference!
# i.e. standard orientation in Gaussian
def ReadXYZ(File):
    f = open(File, "r")
    fc = f.readlines()
    f.close()
    NAtoms = len(fc)
    Coordinates = []
    for line in range(NAtoms):
        Coordinates.append(fc[line].split())
    return Coordinates

def MoleculeFromCoordinates(Atoms, Coordinates):
    """ Takes already read in coordinates and converts
        them into atom objects
    """
    Molecule = dict()
    for Index, Atom in enumerate(Atoms):
        Molecule[Index] = __Atom__(Coordinates[Index,:], Atom)
    return Molecule

def MoleculeFromXYZ(File):
    """ Read in a chemical file format .xyz """
    f = open(File, "r")
    fc = f.readlines()[2:]            # Skip the number of atoms and comment line
    f.close()
    NAtoms = len(fc)
    Molecule = dict()
    for Line in range(NAtoms):
        try:                          # Sometimes the NAtoms is wrong because of a new line at the bottom
            SplitLine = fc[Line].split()
            Symbol = SplitLine[0]                                 # First item is atom symbol
            Coordinates = np.array([SplitLine[1], SplitLine[2], SplitLine[3]])
            Coordinates = Coordinates.astype(np.float)            # convert coordinates to numpy float array
            Molecule[Line] = __Atom__(Coordinates, Symbol)            # Populate dictionary
        except IndexError:
            pass
    return Molecule

def ReadNormalModes(File):
    f = open(File, "r")
    fc = f.readlines()
    f.close()
    NModes = len(fc)
    NormalModes = []
    for line in range(NModes):
        NormalModes.append(fc[line].split())
    return NormalModes

############################### Tools ###############################

# Function that stores a library with all the atomic masses
# and returns it for whatever atom you specify
def Symbol2Mass(Atom):
    # Masses are taken from:
    # http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
    #MassList = {"H": 1.00782503223,
    #            "C": 12.00,
    #            "O": 15.99491461957,
    #            "X": 0.0                                    # Dummy atom
    #            }
#    return MassList[Atom]                                  # Returns mass in amu
    return (masses[Atom] / constants.Avogadro) / 1e3      # Return mass in kg

# Calculates the distance between two atoms
def CalculateDistance(A, B):
    Distance = np.sqrt(np.sum((A.GetCoordinates() - B.GetCoordinates())**2.))
    return Distance

# Uses cartesian points A, B, C to calculate the angle
# formed by A - B - C using vectors. Atom B is the centre of the angle!
def CalculateAngle(A, B, C):
    AB = B.GetCoordinates() - A.GetCoordinates()           # Calculate vector AB
    BC = C.GetCoordinates() - B.GetCoordinates()           # Calculate vector BC
    ABLength = CalculateDistance(A, B)               # Work out magnitude of AB
    BCLength = CalculateDistance(B, C)               # Magnitude of BC
    DotProduct = np.dot(AB, BC)                      # Dot product of AB dot BC
    # Return the angle formed by A - B - C in degrees
    return (np.arccos(DotProduct / (ABLength * BCLength)) * (180. / np.pi))

# Function to return the reduced mass of fragments
# Takes a list of masses in whatever units
def CalculateReducedMass(Masses):
    ReducedMass = 0.
    for index, mass in enumerate(Masses):
        ReducedMass = ReducedMass + (1. / mass)
    return 1. / ReducedMass

############################### Moments of Inertia ###############################

# Here the routines take x,y,z as Angstroms, and mass in kilograms!
# XX diagonal elemnt of I
def IXX(y, z, mass):
    return mass * ((y)**2 + (z)**2)

# ZZ diagonal element of I
def IZZ(x, y, mass):
    return mass * ((x)**2 + (y)**2)

# YY diagonal element of I
def IYY(x, z, mass):
    return mass * ((x)**2 + (z)**2)
def IXY(x, y, mass):
    return mass * (x * y)

def IXZ(x, z, mass):
    return mass * (x * z)

def IYZ(y, z, mass):
    return mass * (y * z)

# Note on this: the standard orientation that Gaussian spits out is already rotated
# to the inertial frame; that means we only need to calculate the diagonal elements
# of the inertia tensor
def OldCalculateInertiaMatrix(Molecule):
    InertiaMatrix = np.zeros((3,3),dtype=float)
    for atom in enumerate(Molecule):
        # Calculate diagonal elements of inertia matrix
        # Enumerate loops through each atom of the molecule;
        # indexes 0:mass, 1:x, 2:y, 3:z
        InertiaMatrix[0,0] = InertiaMatrix[0,0] + IXX(float(atom[1][2]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
        InertiaMatrix[1,1] = InertiaMatrix[1,1] + IYY(float(atom[1][1]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
        InertiaMatrix[2,2] = InertiaMatrix[2,2] + IZZ(float(atom[1][1]), float(atom[1][2]), Symbol2Mass(atom[1][0]))
        # Calculate triangle of off-diagonal elements of inertia matrix
#        InertiaMatrix[0,1] = InertiaMatrix[0,1] + IXY(float(atom[1][1]), float(atom[1][2]), Symbol2Mass(atom[1][0]))
#        InertiaMatrix[0,2] = InertiaMatrix[0,2] + IXZ(float(atom[1][1]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
#        InertiaMatrix[1,2] = InertiaMatrix[1,2] + IYZ(float(atom[1][2]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
        # Symmetrise the matrix
    return (InertiaMatrix + InertiaMatrix.T) / 2.

# Converts the principle moments into rotational constants in 1/cm
# These should agree with the rotational constants in Gaussian
def PMI2ABC(Inertia):
    return constants.h / (8 * (np.pi)**2 * (constants.c * 100.) * Inertia)

############################### Normal mode analysis ###############################

# Take an input cartesian coordinate (x,y,z) as well as the corresponding
# normal mode displacement for that x,y,z
def CalculateDisplacement(x, a):
    return (a * x)**2

# Function to calculate displaced geometries using normal mode displacements
def NormalMode2Cartesian(Molecule, NormalModes, NormalDisplacment):
    DisplacedGeometry = np.zeros((len(Molecule), 3), dtype=float)
    for atom in enumerate(Molecule):
        Displacement = NormalModes[atom[0]] 