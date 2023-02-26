import matplotlib.pyplot as plt
import numpy as np
from Torus import Torus
import csv
import time


T = Torus()

T.loadEigenData("Toric-Codes-Simulator/CSVs/groundEigenvalues.csv",
                "Toric-Codes-Simulator/CSVs/groundEigenstates.csv")


T.printState(T.eigenstates[1])

"""
T.makeHamiltonian()


print("Here we go!")

t0=time.time()

T.makeSparseEigenstatesNaive(k=1)

t1=time.time()

print("Task Successful!")
print("Time elapsed: "+str(t1-t0))

T.makeMWEigendata()


T.saveEigenData("Toric-Codes-Simulator/CSVs/eigenvalues1.csv",
                "Toric-Codes-Simulator/CSVs/eigenstates1.csv")
"""

