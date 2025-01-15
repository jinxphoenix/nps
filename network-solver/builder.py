import numpy as np


class Node:
    _nextID = 0
    def __init__(self):
        self.ID = Node._nextID
        Node._nextID += 1
        self._edges = []
        self.fixedPressure = False 
        self.pressure = 1*10**5 # dmmy variable for development
        self.height  = 0 
        self.demand = 0   
      
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
        
    def setDemand(self,demand):
        self.demand = demand
            

class Edge:
    _nextID = 0
    # will need to pt a lot of the below into a type class (i.e. pipe, pmp, valve etc.)
    def __init__(self,nodeFrom, nodeTo):
        self.ID = Edge._nextID
        self._nodeFrom = nodeFrom
        self._nodeTo = nodeTo
        Edge._nextID += 1 
        self.rougness = 1*10**-5  #dmmy vale for development
        self.length = 100 # dmmy vale for development
        self.hydraulicDiameter = None
        self.area = None
        self._network: Network = None
        self.flowrate = 0.1 # dmmy vale for development
        self.massFlowrate = None
        self.reynolds = None
        self.frictionFactor = 1
        self.resistance = None
        self.gravity = 9.81
        
        self.setDiameter(0.1)  #dmmy vale for development
        
    def attachFromNode(self):
        self._nodeFrom.attachEdge(self)
        
    def attachToNode(self):
        self._nodeTo.attachEdge(self)
        
    def addNetwork(self,network):
        self._network = network
        
    def setDiameter(self,diameter):
        self.hydraulicDiameter = diameter
        self.calculateArea()
        
        
    def calculateArea(self):
        self.area = np.pi*self.hydraulicDiameter/4
        
    def updateReynolds(self):
        #self.reynolds = self.massFlowrate * self.hydraulicDiameter /(self._network._fluid.viscosity * self.area )
        
        self.reynolds = self.flowrate * self.hydraulicDiameter * self._network._fluid.density /(self._network._fluid.viscosity * self.area)
         
    def updatefrictionFactor(self):
        self.updateReynolds()
        if self.reynolds > 2000:
            a = self.rougness/(3.7*self.hydraulicDiameter)
            b = 2.51 / (self.rougness*np.sqrt(self.frictionFactor))
            c = -2*np.log10(a + b)
            self.frictionFactor = (1/c)**2
        else:
            self.frictionFactor = 64/self.reynolds
            
    def updateResistance(self):
        self.updatefrictionFactor()
        self.resistance = (8 * self.frictionFactor * self.length) / (np.pi**2 * self.gravity *  self.hydraulicDiameter**5)
        return self.resistance            

class Network:
    def __init__(self):
        self._nodes: list[Node] = []
        self._edges: list[Edge]  = []
        self._unodes: list[Node] = []
        self._fnodes: list[Node] = []
        self._fluid: Fluid = None
                
    def addEdge(self,edge):
        if edge not in self._edges:
            self._edges.append(edge)
            edge.addNetwork(self)
            
    def addNode(self,node):
        if node not in self._nodes:
            self._nodes.append(node)
            if node.fixedPressure == False:
                self._unodes.append(node)
            else:
                self._fnodes.append(node)
    
    def addFluid(self,fluid):
        self._fluid = fluid
                    
    def generateConnectivityMatrix(self):
        # this is the A1 amtri
        self.connectivity = np.zeros((len(self._unodes),len(self._edges)))
        i = 0
        for n in self._unodes:
            j = 0
            for e in self._edges:
                if e._nodeFrom == n:
                    self.connectivity[i, j] = 1
                    
                elif e._nodeTo == n:
                    self.connectivity[i, j] = -1
                j += 1
            i += 1
        self.connectivity = np.transpose(self.connectivity)    
                
    def generateFixedMatrix(self):
        # this is the A2 matri
        self.fixed = np.zeros((len(self._fnodes),len(self._edges)))
        i = 0
        for n in self._fnodes:
            j = 0
            for e in self._edges:
                if e._nodeFrom == n:
                    self.fixed[i, j] = 1
                    
                elif e._nodeTo == n:
                    self.fixed[i, j] = -1   
                j += 1
            i += 1
        self.fixed = np.transpose(self.fixed)
        
    def generateResistanceMatrix(self):
        self.resistances = np.zeros((len(self._edges),len(self._edges)))
        i = 0
        for e in self._edges:
            self.resistances[i,i] = e.updateResistance()
            i += 1
             
    def createFlowVector(self):
        self.flowVector = np.zeros((len(self._edges),1))
        i = 0
        for e in self._edges:
            self.flowVector[i] = e.flowrate
            i += 1
            
    def createDemandVector(self):
        self.demandVector = np.zeros((len(self._edges),1))
        i = 0
        for n in self._unodes:
            self.demandVector[i] = n.demand
            i += 1
    
    def generateGMatrix(self):
        # this is the G matrix
        self.G = np.multiply(self.resistances,np.abs(self.flowVector[:]))
                     

class Fluid:
    def __init__(self,density,viscosity):
        self.density = density
        self.viscosity = viscosity
        

net0 = Network()     

fl0 = Fluid(1000,1*10**-3) 

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()

e01 = Edge(n0,n1)
e12 = Edge(n1,n2)
e13 = Edge(n1,n3)

n2.setDemand(10/3600)
n3.setDemand(5/3600)
n0.setFixedPressure(10)

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

net0.generateConnectivityMatrix()
print(net0.connectivity)

net0.generateFixedMatrix()
print(net0.fixed)

net0.generateResistanceMatrix()
print(net0.resistances)

net0.createFlowVector()
print(net0.flowVector)

net0.generateGMatrix()
print(net0.G)

net0.createDemandVector()
print(net0.demandVector)