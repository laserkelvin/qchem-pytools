# qchem-pytools
Python routines to turn quantum/physical chemistry into something meaningful.

The code is written in Python 3.5, and should be backwards (2.7) compatible.

Finally got around to writing some PIP ready code.

## Installation

1. Clone this git
2. Extract the files somewhere
3. Navigate to the root folder containing setup.py
4. Run `pip install .`

---

# Routine descriptions

## cfour_parse.py

This is a parser I wrote to take the most important output from CFOUR calculations. All of the information is extracted and organised into a dictionary, which can be dumped as JSON, HTML, or just used as a Python dictionary. 

This is written with two things in mind: __reproducibility__ and __time saving__. I want to be able to look at this report and figure out what was done and the results are _quickly_.

The routines are kept in a class; `OutputFile`.

The dictionary is organised into properties and information about a calculation. The method `parse` will open a target file, and read in line-by-line (to minimise memory usage). The way this parser works is simple string matching, and flagging.

The flagging is done because I'm reading in line-by-line, rather than keeping the whole file in memory which would allow indexing. Since the file is streamed, I have no way of keeping track of line numbers, which would be useful when pulling multi-line information. Instead, Booleans are used to keep track of when the parser records information, and when to stop.

An illustration with extracting rotational constants: the rotational constants are kept in the line __after__ the string "Rotational constants". When this line is matched, the `RotFlag` variable is used to track when to read in constants.

1. Current line contains "Rotational constants (in MHz)"
2. `RotFlag` is set to `True`, and we move on.
3. Every time a line is read, RotFlag is checked:
    1. If `True`, we read in the rotational constants. `RotFlag` is set to `False`.
    2. If `False`, check the other `if` cases.

---

## figure_settings.py

These are routines I wrote that will let me quickly plot something up.
