import matplotlib.pyplot as plt
import numpy as np
from Torus import Torus

S = Torus(subdivisions=2)

S.makeHamiltonian()

S.makeEigenstates()

dataMinus8=[S.eigenstates[i] for i in range(256) if S.eigenvalues[i]==-8]
dataMinus4=[S.eigenstates[i] for i in range(256) if S.eigenvalues[i]==-4]
data0=[S.eigenstates[i] for i in range(256) if S.eigenvalues[i]==0]
data4=[S.eigenstates[i] for i in range(256) if S.eigenvalues[i]==4]
data8=[S.eigenstates[i] for i in range(256) if S.eigenvalues[i]==8]

valsMinus8=[S.compute_Q_ptrace(data) for data in dataMinus8]
valsMinus4=[S.compute_Q_ptrace(data) for data in dataMinus4]
vals0=[S.compute_Q_ptrace(data) for data in data0]
vals4=[S.compute_Q_ptrace(data) for data in data4]
vals8=[S.compute_Q_ptrace(data) for data in data8]

plt.plot(vals8, np.zeros_like(vals8)+8, 'x',label="Eigenvalue 8")
plt.plot(vals4, np.zeros_like(vals4)+4, 'x',label="Eigenvalue 4")
plt.plot(vals0, np.zeros_like(vals0), 'x',label="Eigenvalue 0")
plt.plot(valsMinus4, np.zeros_like(valsMinus4)-4, 'x',label="Eigenvalue -4")
plt.plot(valsMinus8, np.zeros_like(valsMinus8)-8, 'x',label="Eigenvalue -8")
plt.legend()
plt.show()

print(len(dataMinus8))

print(len(dataMinus4))

print(len(data0))

print(len(data4))

print(len(data8))