from setuptools import setup

setup(
    name="qchem_pytools",
    version="0.3",
    description="A set of Python routines for quantum/physical chemistry analysis",
    author="Kelvin Lee",
    packages=["qchem_pytools"],
    author_email="kin.long.kelvin.lee@gmail.com",
    install_requires=[
            "numpy",
            "plotly",
            "pandas",
            "scipy",
            "peakutils",
            "seaborn",
            #"xhtml2pdf",
            "colorlover",
            "matplotlib"
    ]
)
