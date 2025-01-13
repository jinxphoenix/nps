import numpy as np


class Node:
    _nextID = 0
    def __init__(self):
        self.ID = Node._nextID
        Node._nextID += 1
        self._edges = []
        self.fixedPressure = False 
        self.pressure = 0
        self.height  = 0    
      
    def attachEdge(self,edge):
        if edge not in self._edges:
            self._edges.append(edge)
    
    def detachEdge(self,edge):
        try:          
            self._edges.remove(edge)
            
        except ValueError:
            pass
    
    def setFixedPressure(self,pressure):
        self.fixedPressure = True
        self.pressure = pressure
        
    def setHeight(self,height):
        self.height = height

class Edge:
    _nextID = 0
    def __init__(self,nodeFrom, nodeTo):
        self.ID = Edge._nextID
        self._nodeFrom = nodeFrom
        self._nodeTo = nodeTo
        Edge._nextID += 1 
        
    def attachFromNode(self):
        self._nodeFrom.attachEdge(self)
        
    def attachToNode(self):
        self._nodeTo.attachEdge(self)
                    

class Network:
    def __init__(self):
        self._nodes = []
        self._edges = []
        self._unodes = []
        self._fnodes = []
        self._fluid = None
                
    def addEdge(self,edge):
        if edge not in self._edges:
            self._edges.append(edge)
            
    def addNode(self,node):
        if node not in self._nodes:
            self._nodes.append(node)
            if node.fixedPressure == False:
                self._unodes.append(node)
            else:
                self._fnodes.append(node)
    
    def addFluid(self,fluid):
        self._fluid = fluid
                    
    def generateIncidenceMatrix(self):
        self.incidence = np.zeros((len(self._nodes),len(self._edges)))
        for n in self._nodes:
            for e in self._edges:
                if e._nodeFrom == n:
                    self.incidence[n.ID, e.ID] = -1
                    
                elif e._nodeTo == n:
                    self.incidence[n.ID, e.ID] = 1
                      
    
    def generateAdjacencyMatrix(self):
        self.adjacency = np.zeros((len(self._nodes),len(self._nodes)))
        for n in self._nodes:
            for e in self._edges:
                if e._nodeFrom == n:
                    self.adjacency[n.ID, e._nodeTo.ID] = 1
                    
                if e._nodeTo == n:
                    self.adjacency[n.ID, e._nodeFrom.ID] = 1

    def generateConnectivityMatrix(self):
        self.connectivity = np.zeros((len(self._unodes),len(self._edges)))
        i = 0
        for n in self._unodes:
            j = 0
            for e in self._edges:
                if e._nodeFrom == n:
                    self.connectivity[n.ID, e.ID] = -1
                    
                elif e._nodeTo == n:
                    self.connectivity[n.ID, e.ID] = 1
                j += 1
            i += 1
                
    def generateFixedMatrix(self):
        self.fixed = np.zeros((len(self._fnodes),len(self._edges)))
        i = 0
        for n in self._fnodes:
            j = 0
            for e in self._edges:
                if e._nodeFrom == n:
                    self.fixed[i, j] = -1
                    
                elif e._nodeTo == n:
                    self.fixed[i, j] = 1   
                j += 1
            i += 1                    
                     

class Fliud:
    def __init__(self,density,viscosity):
        self.density = density
        self.viscosity = viscosity
        

net0 = Network()     

fl0 = Fliud(1000,1*10**-3) 

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()

e01 = Edge(n0,n1)
e12 = Edge(n1,n2)
e13 = Edge(n1,n3)

n2.setFixedPressure(1)
n3.setFixedPressure(1)

net0.addFluid(fl0)
net0.addNode(n0)
net0.addNode(n1)
net0.addNode(n2)
net0.addNode(n3)

net0.addEdge(e01)
net0.addEdge(e12)
net0.addEdge(e13)

print(len(net0._edges))
print(len(net0._nodes))
print(len(net0._unodes))
print(len(net0._fnodes))

net0.generateIncidenceMatrix()
print(net0.incidence)

net0.generateAdjacencyMatrix()
print(net0.adjacency)

net0.generateConnectivityMatrix()
print(net0.connectivity)

net0.generateFixedMatrix()
print(net0.fixed)