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
                indicies+="Z"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Get the matrix for applying Pauli(Y) to neighbors of a face
    def B(self,face):

        indicies=""

        for edge in self.edges:
            if edge in self.faceNeighbors(face):
                indicies+="X"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Get the operator corresponding to a homology class
    def homologyOperator(self,homology):

        indicies=""

        if homology=="id":
            Xs=[]

        if homology=="alpha":
            Xs=[self.makeHash([(2,0),(2,1)]),
                self.makeHash([(2,1),(2,2)]),
                self.makeHash([(2,2),(2,0)])]
        
        if homology=="beta":
            Xs=[self.makeHash([(0,2),(1,2)]),
                self.makeHash([(1,2),(2,2)]),
                self.makeHash([(2,2),(0,2)])]

        if homology=="alphabeta":
            Xs=[self.makeHash([(2,0),(2,1)]),
                self.makeHash([(2,1),(2,2)]),
                self.makeHash([(2,2),(2,0)]),
                self.makeHash([(0,2),(1,2)]),
                self.makeHash([(1,2),(2,2)]),
                self.makeHash([(2,2),(0,2)])]
            
        for edge in self.edges:
            if edge in Xs:
                indicies+="X"
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
    def makeEigenstatesNaive(self):

        mat=self.hamiltonian.to_matrix()

        data=np.linalg.eig(mat)

        self.eigenvalues=np.round(data[0])
        self.eigenstates=np.transpose(data[1])
  
    #Get the first k eigenstates and eigenvalues, stored in self.eigenvalues and self.eigenstates
    def makeSparseEigenstatesNaive(self,k):

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
    
    #Print a (pure) state vector in a reasonable format
    def printPureState(self,ket):
        pretty=""
        pretty+=" __A__  __B__ __C__"+"\n"
        pretty+="|     |      |      |"+"\n"
        pretty+="D     E      F      D"+"\n"
        pretty+="|     |      |      |"+"\n"
        pretty+=" __G__ __H__   __I__"+"\n"
        pretty+="|     |      |      |"+"\n"
        pretty+="J     K      L      J"+"\n"
        pretty+="|     |      |      |"+"\n"
        pretty+=" __M__ __N__   __O__"+"\n"
        pretty+="|     |      |      |"+"\n"
        pretty+="P     Q      R      P"+"\n"
        pretty+="|     |      |      |"+"\n"
        pretty+=" __A__ __B__   __C__"

        pretty=pretty.replace("P",str(ket[0]))
        pretty=pretty.replace("Q",str(ket[1]))
        pretty=pretty.replace("R",str(ket[2]))
        pretty=pretty.replace("J",str(ket[3]))
        pretty=pretty.replace("K",str(ket[4]))
        pretty=pretty.replace("L",str(ket[5]))
        pretty=pretty.replace("D",str(ket[6]))
        pretty=pretty.replace("E",str(ket[7]))
        pretty=pretty.replace("F",str(ket[8]))
        pretty=pretty.replace("A",str(ket[9]))
        pretty=pretty.replace("B",str(ket[10]))
        pretty=pretty.replace("C",str(ket[11]))
        pretty=pretty.replace("M",str(ket[12]))
        pretty=pretty.replace("N",str(ket[13]))
        pretty=pretty.replace("O",str(ket[14]))
        pretty=pretty.replace("G",str(ket[15]))
        pretty=pretty.replace("H",str(ket[16]))
        pretty=pretty.replace("I",str(ket[17]))

        print(pretty)

    #Present a state vector in a reasonable format
    def printState(self,ket):
        L=int(np.log2(len(ket)))
        readable={}

        for pureState in range(len(ket)):
            if ket[pureState]!=0:
                zeros = "0"*(L-len(bin(pureState)[2:]))
                readable[zeros+bin(pureState)[2:]]=ket[pureState]

        for item in readable:
            self.printPureState(item)
            print(readable[item])

    #Get the values and multiplicities of a vector
    def getVals(self,ket):
        values={}
        for val in ket:
            if val not in values:
                values[val]=1
            else:
                values[val]+=1
        return values
    
    #Plot data
    def plotEigendata(self):
    
        for eigenval in sorted(list(set(self.eigenvalues)),reverse=True):
            data=[self.eigenMW[i] for i in range(len(self.eigenvalues)) if self.eigenvalues[i]==eigenval]
            print(str(eigenval)+": "+str(len(data)))
            plt.plot(data, np.zeros_like(data)+eigenval, 'x',label="Eigenvalue "+str(int(eigenval)))
        
        plt.legend()
        plt.show()

    #Save data
    def saveEigenData(self, eigenvalueCSV, eigenstateCSV):

        with open(eigenvalueCSV,"w") as f:

            writer = csv.writer(f)

            for eigenvalue in self.eigenvalues:
                writer.writerow([eigenvalue])

        with open(eigenstateCSV,"w") as f:

            writer = csv.writer(f)

            for eigenstate in self.eigenstates:
                writer.writerow(eigenstate)
        
    #Load data
    def loadEigenData(self, eigenvalueCSV,eigenstateCSV):

        with open(eigenvalueCSV,"r") as f:

            reader = csv.reader(f)

            for eigenvalue in reader:
                self.eigenvalues+=[int(np.real(complex(eigenvalue[0])))]

        with open(eigenstateCSV,"r") as f:

            reader = csv.reader(f)

            for eigenstate in reader:
                self.eigenstates+=[[complex(eigenstate[i]) for i in range(len(eigenstate))]]
        

