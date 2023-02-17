#Imports
import numpy as np
import qutip

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister,BasicAer
from qiskit.compiler import transpile
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.quantum_info import process_fidelity
from qiskit.extensions import RXGate, XGate, CXGate

#Main "Torus" class
class Torus:
    #Initialize the Torus
    def __init__(self,subdivisions=3):
        
        self.N = subdivisions

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

        matrix=np.array([1])

        for edge in self.edges:
            if edge in self.vertexNeighbors(vertex):
                matrix=np.kron(matrix,Pauli("X"))
            else:
                matrix=np.kron(matrix,Pauli("I"))
        
        return matrix

    #Get the matrix for applying Pauli(Y) to neighbors of a face
    def B(self,face):

        matrix=np.array([1])

        for edge in self.edges:
            if edge in self.faceNeighbors(face):
                matrix=np.kron(matrix,Pauli("Y"))
            else:
                matrix=np.kron(matrix,Pauli("I"))
        
        return matrix

    #Make the hamiltonian of the torus, stored in self.hamiltonian
    def makeHamiltonian(self):

        self.hamiltonian=np.zeros([2**(2*self.N**2),2**(2*self.N**2)])

        for vertex in self.verticies:
            self.hamiltonian=self.hamiltonian-self.A(vertex)
        
        for face in self.faces:
            self.hamiltonian=self.hamiltonian-self.B(face)
        
    #Get the eigenstates and eigenvalues, stored in self.eigenvalues and self.eigenstates
    def makeEigenstates(self):

        data=np.linalg.eig(self.hamiltonian)

        self.eigenvalues=np.round(data[0])
        self.eigenstates=data[1]

    #Compute the Mayer-Wallash entropy of a state ket
    def compute_Q_ptrace(self,ket):
        N=int(np.log2(len(ket)))
        ket = qutip.Qobj(ket, dims=[[2]*(N), [1]*(N)]).unit()
        entanglement_sum = 0
        for k in range(N):
            rho_k_sq = ket.ptrace([k])**2
            entanglement_sum += rho_k_sq.tr()  
    
        Q = 2*(1 - (1/N)*entanglement_sum)
        return Q

