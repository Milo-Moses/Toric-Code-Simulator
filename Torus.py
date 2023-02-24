#Imports
import numpy as np
import qutip
import csv


from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister,BasicAer
from qiskit.compiler import transpile
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.quantum_info import process_fidelity
from qiskit.extensions import RXGate, XGate, CXGate
from qiskit.opflow.primitive_ops import PauliSumOp
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt


#Main "Torus" class
class Torus:
    #Initialize the Torus
    def __init__(self,subdivisions=3,eigenstates=[],eigenvalues=[],eigenMW=[]):
        
        self.N = subdivisions

        self.eigenstates=eigenstates

        self.eigenvalues=eigenvalues

        self.eigenMW=eigenMW

        self.verticies = [
            self.makeHash([(i,j)])
            for i in range(self.N)
            for j in range(self.N)
        ]

        self.edges = [
            self.makeHash([(i,j),(i+1,j)])
            for i in range(self.N)
            for j in range(self.N)
        ]+[
            self.makeHash([(i,j),(i,j+1)])
            for i in range(self.N)
            for j in range(self.N)
        ]

        self.edgesDict = {k: edge for edge, k in enumerate(self.edges)}

        self.faces = [
            self.makeHash([(i,j),(i+1,j),(i,j+1),(i+1,j+1)])
            for i in range(self.N)
            for j in range(self.N)
        ]

        #The internal state of our Torus, index by edges via edgesDict
        self.state = QuantumCircuit(0)

        for edge in self.edges:
            self.state.add_register(QuantumRegister(1,edge))

    #Get the neighboring edges of a vertex        
    def vertexNeighbors(self,vertex):
        neighbors=[]

        for edge in self.edges:
            if vertex in edge.split(" "):
                neighbors+=[edge]
        
        return neighbors

    #Get the neighboring edges of a face
    def faceNeighbors(self,face):
        neighbors=[]

        for edge in self.edges:
            
            if (edge.split(" ")[0] in face.split(" ")) and (edge.split(" ")[1] in face.split(" ")):
                neighbors+=[edge]
        
        return neighbors
    
    #Hash a tuple (i,j) into a hashable item for edgesDict
    def makeHash(self,tuples):
        hashed=""
    
        for pair in tuples:
            hashed+=str(pair[0]%self.N)
            hashed+="."
            hashed+=str(pair[1]%self.N)
            hashed+=" "
        
        hashed=hashed[:-1]

        return hashed

    #Get the matrix for applying Pauli(X) to neighbors of a vertex         
    def A(self,vertex):

        indicies=""

        for edge in self.edges:
            if edge in self.vertexNeighbors(vertex):
                indicies+="X"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Get the matrix for applying Pauli(Y) to neighbors of a face
    def B(self,face):

        indicies=""

        for edge in self.edges:
            if edge in self.faceNeighbors(face):
                indicies+="Y"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Make the hamiltonian of the torus, stored in self.hamiltonian
    def makeHamiltonian(self):

        self.hamiltonian=PauliSumOp(SparsePauliOp([2*(self.N**2)*"I"],[0]))

        for vertex in self.verticies:
            self.hamiltonian=self.hamiltonian-self.A(vertex)
        
        for face in self.faces:
            self.hamiltonian=self.hamiltonian-self.B(face)
    
    #Get the eigenstates and eigenvalues, stored in self.eigenvalues and self.eigenstates
    def makeEigenstates(self):

        mat=self.hamiltonian.to_matrix()

        data=np.linalg.eig(mat)

        self.eigenvalues=np.round(data[0])
        self.eigenstates=np.transpose(data[1])
  
    #Get the first k eigenstates and eigenvalues, stored in self.eigenvalues and self.eigenstates
    def makeSparseEigenstates(self,k):

        mat=self.hamiltonian.to_spmatrix()

        data=eigs(mat,k=k)

        self.eigenvalues=np.round(data[0])
        self.eigenstates=np.transpose(data[1])

    #Compute the Mayer-Wallash entropy of a state ket
    def MW(self,ket):
        N=int(np.log2(len(ket)))
        ket = qutip.Qobj(ket, dims=[[2]*(N), [1]*(N)]).unit()
        entanglement_sum = 0
        for k in range(N):
            rho_k_sq = ket.ptrace([k])**2
            entanglement_sum += rho_k_sq.tr()  
    
        Q = 2*(1 - (1/N)*entanglement_sum)
        return Q

    #Compute the Mayer-Wallash entropy of all the eigenstates
    def makeMWEigendata(self):
        self.eigenMW=[]

        for eigenstate in self.eigenstates:
            self.eigenMW+=[self.MW(eigenstate)]
    
    #Plot data
    def plotEigendata(self):
    
        for eigenval in sorted(list(set(self.eigenvalues)),reverse=True):
            data=[self.eigenMW[i] for i in range(len(self.eigenvalues)) if self.eigenvalues[i]==eigenval]

            plt.plot(data, np.zeros_like(data)+eigenval, 'x',label="Eigenvalue "+str(int(eigenval)))
        
        plt.legend()
        plt.show()

    #Save data
    def saveEigenData(self, eigenvalueCSV, eigenMWCSV, eigenstateCSV):

        with open(eigenvalueCSV,"w") as f:

            writer = csv.writer(f)

            for eigenvalue in self.eigenvalues:
                writer.writerow([eigenvalue])

        with open(eigenMWCSV,"w") as f:

            writer = csv.writer(f)

            for MW in self.eigenMW:
                writer.writerow([MW])

        with open(eigenstateCSV,"w") as f:

            writer = csv.writer(f)

            for eigenstate in self.eigenstates:
                writer.writerow(eigenstate)
        
    #Load data
    def loadEigenData(self, eigenvalueCSV, eigenMWCSV,eigenstateCSV):

        with open(eigenvalueCSV,"r") as f:

            reader = csv.reader(f)

            for eigenvalue in reader:
                self.eigenvalues+=[int(np.real(complex(eigenvalue[0])))]

        with open(eigenMWCSV,"r") as f:

            reader = csv.reader(f)

            for MW in reader:
                self.eigenMW+=[float(MW[0])]

        with open(eigenstateCSV,"r") as f:

            reader = csv.reader(f)

            for eigenstate in reader:
                self.eigenMW+=[complex(eigenstate[i]) for i in range(len(eigenstate))]
        

