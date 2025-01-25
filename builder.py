
from Pipe import Pipe
import numpy as np


class Node:
    _nextID = 0
    def __init__(self):
        self.ID = Node._nextID
        Node._nextID += 1
        self._edges = []
        
        self.pressure = 1 # dmmy variable for development
        self.previousPressure = self.pressure*0.9 # delib off set for intial itteration.
        
        self.height  = 0 
        self.demand = 0   
        self.fixedPressure = False 
        
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
    
    def __init__(self,nodeFrom, nodeTo, model):
        self.ID = Edge._nextID
        self._nodeFrom = nodeFrom
        self._nodeTo = nodeTo
        Edge._nextID += 1 
        self._network: Network = None
        
        self.flowrate = 0.1 # dmmy vale for development
        self.previousFlowrate = self.flowrate - 1 # Initiallised to a value of not the prev flow rate so first itteration works
        self.massFlowrate = None
        
        self.model: Pipe = model
                
        
    def attachFromNode(self):
        self._nodeFrom.attachEdge(self)
        
    def attachToNode(self):
        self._nodeTo.attachEdge(self)
        
    def addNetwork(self,network):
        self._network = network       

    def updateResistance(self):
        return self.model.updateResistance(self.flowrate)


class Network:
    def __init__(self):
        self._nodes: list[Node] = []
        self._edges: list[Edge]  = []
        self._unodes: list[Node] = []
        self._fnodes: list[Node] = []
        self._fluid: Fluid = None
        self.itterationCount = 0
        self.flowTolerance = 1e-4
        self.pressureTolerance = 1e-4
                
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
        # this is the A1 matrix num edges x bum unkown pressure nodes
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
    
    def createUnkownHeadVector(self):
        self.unkownHeadVector = np.zeros((len(self._unodes),1))
        i = 0
        for n in self._unodes:
            self.unkownHeadVector[i] = n.pressure
            i += 1
             
    def createFlowVector(self):
        self.flowVector = np.zeros((len(self._edges),1))
        i = 0
        for e in self._edges:
            self.flowVector[i] = e.flowrate
            i += 1
            
    def createDemandVector(self):
        self.demandVector = np.zeros((len(self._unodes),1))
        i = 0
        for n in self._unodes:
            self.demandVector[i] = n.demand
            i += 1
    
    def generateGMatrix(self):
        # this is the G matrix
        
        # regenerate the resistance
        self.generateResistanceMatrix()
        
        # then create the G matrix
        self.G = np.multiply(self.resistances,np.abs(self.flowVector[:]))
        
    def generateNMatrix(self):
        """ generates an N matrix based on the exponents required in the Jacobian for each edge
            For each edges, it examines what exponent is used for flow (darcey = 2, hazen = 1.852)
            and creates a diagonal matrix for each edge with that value of N.
            
            This is used as part of the jacobian (first derivative) and is thus used by 
        """
        self.N = np.eye(len(self._edges),len(self._edges))
        
    def createFixedHeadVector(self):
        self.fixedHeadVector = np.zeros((len(self._fnodes),1))
        i = 0
        for n in self._fnodes:
            self.fixedHeadVector[i] = n.pressure
            i += 1
    
    def calculateFlowResidual(self):
        self.flowResidual = np.zeros((len(self._edges),1))
        i = 0
        for e in self._edges:
            self.flowResidual[i] = e.flowrate - e.previousFlowrate
            i += 1
        
        max_res = np.abs(np.max(self.flowResidual))
        print(f'max flow residual is: {max_res}')
        return max_res
    
    def calculatePressureResidual(self):
        self.pressureResidual = np.zeros((len(self._unodes),1))
        i = 0
        for n in self._unodes:
            self.pressureResidual[i] = n.pressure - n.previousPressure
            i += 1
        
        max_res = np.abs(np.max(self.pressureResidual))
        print(f'max pressure residual is: {max_res}')
        return max_res
        
        
                
    def initialise(self):
        """Initialises all the matrices and vectors required for the calculation.
        """
        
        self.generateConnectivityMatrix()
        self.generateFixedMatrix()
        self.createFlowVector()
        self.createDemandVector()
        self.generateNMatrix()
        self.createUnkownHeadVector()
        self.createFixedHeadVector()
    
    def itterateHead(self):
        # itterate the head value
        
        # first store the previous itterations pressures
        i = 0
        for n in self._unodes:
            n.previousPressure = n.pressure
            i += 1     
        
        # note n = 2 at the moment. This needs to be refactored to use the n matrix and N matrix function.
        schur = np.linalg.inv(self.connectivity.transpose() @ np.linalg.inv(self.G) @ self.connectivity)
        other1 = self.connectivity.transpose() @ ((1-2) * self.flowVector - np.linalg.inv(self.G) @ (self.fixed @ self.fixedHeadVector))
        other2 =  - 2 * self.demandVector
        other = (other1 + other2)
        self.unkownHeadVector = schur @ other
        print(self.unkownHeadVector)
        
        i = 0
        for n in self._unodes:
            n.pressure = self.unkownHeadVector[i]
            i += 1
        
    def itterateFlow(self):
        # itterate flow (after itterating head)
        
        i = 0
        for e in self._edges:
            e.previousFlowrate = e.flowrate
            i += 1
        
        # note n = 2 at the moment. This needs to be refactored to use the n matrix and N matrix function.                   
        self.flowVector = (1-1/2)*self.flowVector + (1/2)*np.linalg.inv(self.G) @ (self.connectivity @ self.unkownHeadVector + self.fixed @ self.fixedHeadVector)
        print(self.flowVector)
        
        i = 0
        for e in self._edges:
            e.flowrate = self.flowVector[i]
            i += 1
        
    def itterate(self,ittLimit):
        itterations = 0
        while ( itterations<ittLimit and ((self.calculateFlowResidual() > self.flowTolerance) or (self.calculatePressureResidual() > self.pressureTolerance))):
            print(f'Itterations: {itterations}')
            net0.generateGMatrix()       
            self.itterateHead()
            self.itterateFlow()
            itterations += 1
            print(net0.flowVector)
            print(net0.unkownHeadVector)    
            print(net0.fixedHeadVector)
        
            

class Fluid:
    def __init__(self,density,viscosity):
        self.density = density
        self.viscosity = viscosity

# Set up the network
net0 = Network()     

fl0 = Fluid(1000,1*10**-3) 

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()
n4 = Node()

p = Pipe(1*10**-5,100,0.05)

e01 = Edge(n0,n1,p)
e12 = Edge(n1,n2,p)
e13 = Edge(n1,n3,p)
e14 = Edge(n1,n4,p)

n2.setDemand(10/3600)
n3.setDemand(5/3600)
# n4.setFixedPressure(10)
n4.setDemand(5/3600)
n0.setFixedPressure(20)

net0.addFluid(fl0)
net0.addNode(n0)
net0.addNode(n1)
net0.addNode(n2)
net0.addNode(n3)
net0.addNode(n4)

net0.addEdge(e01)
net0.addEdge(e12)
net0.addEdge(e13)
net0.addEdge(e14)

# initialise and run the simulation

net0.initialise()
net0.itterate(10)