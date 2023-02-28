#Imports
import numpy as np
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

        self.makeHamiltonian()

        self.loadEigenData("CSVs/groundEigenstates.csv")

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
    def makeHash(self,tuples,hashed=False):

        reducedTuple=[(pair[0]%self.N,pair[1]%self.N) for pair in tuples]

        reducedTuple=sorted(reducedTuple)

        hashed=""
    
        for pair in reducedTuple:
            hashed+=str(pair[0])
            hashed+="."
            hashed+=str(pair[1])
            hashed+=" "
        
        hashed=hashed[:-1]

        return hashed

    #Apply Pauli(X) to an edge      
    def X(self,edge0,hashed=False):

        if not hashed:
            edge0=self.makeHash(edge0)

        indicies=""

        for edge in self.edges:
            if edge==edge0:
                indicies+="X"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Apply Pauli(Z) to an edge      
    def Z(self,edge0,hashed=False):

        if not hashed:
            edge0=self.makeHash(edge0)

        indicies=""

        for edge in self.edges:
            if edge==edge0:
                indicies+="Z"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))
    
    #Get the matrix for applying Pauli(Z) to neighbors of a vertex         
    def A(self,vertex,hashed=False):

        if not hashed:
            vertex=self.makeHash(vertex)

        indicies=""

        for edge in self.edges:
            if edge in self.vertexNeighbors(vertex):
                indicies+="Z"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Get the matrix for applying Pauli(X) to neighbors of a face
    def B(self,face,hashed=False):

        if not hashed:
            face=self.makeHash(face)

        indicies=""

        for edge in self.edges:
            if edge in self.faceNeighbors(face):
                indicies+="X"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))

    #Get the operator corresponding to a homology class
    def homologyXOperator(self,homology):

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

    #Get the operator corresponding to a homology class
    def homologyZOperator(self,homology):

        indicies=""

        if homology=="id":
            Zs=[]

        if homology=="alpha":
            Zs=[self.makeHash([(1,0),(2,0)]),
                self.makeHash([(1,1),(2,1)]),
                self.makeHash([(1,2),(2,2)])]
        
        if homology=="beta":
            Zs=[self.makeHash([(0,1),(0,2)]),
                self.makeHash([(1,1),(1,2)]),
                self.makeHash([(2,1),(2,2)])]

        if homology=="alphabeta":
            Zs=[self.makeHash([(1,0),(2,0)]),
                self.makeHash([(1,1),(2,1)]),
                self.makeHash([(1,2),(2,2)]),
                self.makeHash([(0,1),(0,2)]),
                self.makeHash([(1,1),(1,2)]),
                self.makeHash([(2,1),(2,2)])]
            
        for edge in self.edges:
            if edge in Zs:
                indicies+="Z"
            else:
                indicies+="I"
        
        return PauliSumOp(SparsePauliOp([indicies],[1]))
    
    #Make the hamiltonian of the torus, stored in self.hamiltonian
    def makeHamiltonian(self):

        self.hamiltonian=PauliSumOp(SparsePauliOp([2*(self.N**2)*"I"],[0]))

        for vertex in self.verticies:
            self.hamiltonian=self.hamiltonian-self.A(vertex,hashed=True)
        
        for face in self.faces:
            self.hamiltonian=self.hamiltonian-self.B(face,hashed=True)
    
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
      
    #Load data
    def loadEigenData(self,eigenstateCSV):

        eigenstates=[]

        with open(eigenstateCSV,"r") as f:

            reader = csv.reader(f)

            for eigenstate in reader:
                eigenstates+=[[complex(eigenstate[i]) for i in range(len(eigenstate))]]

        self.ground={}
        self.ground["0"]=eigenstates[0]
        self.ground["alpha"]=eigenstates[1]
        self.ground["beta"]=eigenstates[2]
        self.ground["alphabeta"]=eigenstates[3]

        

        

