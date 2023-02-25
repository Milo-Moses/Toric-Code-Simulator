import matplotlib.pyplot as plt
import numpy as np
from Torus import Torus
import csv
import time


T = Torus()


T.makeHamiltonian()


print("Here we go!")

t0=time.time()

T.makeSparseEigenstatesNaive(k=1)

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

values={}

for item in T.eigenstates:

    if item in values:
        values[item]+=1
    
    else:
        values[item]=1

for item in values:
    print(str(item)+": "+str(values[item]))

"""