import matplotlib.pyplot as plt
import numpy as np
from Torus import Torus
import csv
import time


T = Torus(subdivisions=3)

T.makeHamiltonian()

print("Here we go!")

t0=time.time()

T.makeSparseEigenstates(k=3)

t1=time.time()

print("Task Successful!")
print("Time elapsed: "+str(t1-t0))

T.makeMWEigendata()

T.plotEigendata()

T.saveEigenData("Surface_Codes/CSVs/eigenvalues1.csv",
                "Surface_Codes/CSVs/eigenMW1.csv",
                "Surface_Codes/CSVs/eigenstates1.csv")

"""
T.loadEigenData("Surface_Codes/CSVs/eigenvalues1.csv",
                "Surface_Codes/CSVs/eigenMW1.csv",
                "Surface_Codes/CSVs/eigenstates1.csv")

T.plotEigendata()
"""