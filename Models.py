
from enum import Enum
import numpy as np

from test import N

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
    
    def __init__(self,nodeFrom, nodeTo):
        self.ID = Edge._nextID
        self._nodeFrom: Node = nodeFrom
        self._nodeTo: Node = nodeTo
        Edge._nextID += 1 
        self._network: Network = None
        
        self.status = None
        self.flowrate = 0.1 # dmmy vale for development
        self.previousFlowrate = self.flowrate - 1 # Initiallised to a value of not the prev flow rate so first itteration works
        self.N_value = 1
                  
    def attachFromNode(self):
        self._nodeFrom.attachEdge(self)
        
    def attachToNode(self):
        self._nodeTo.attachEdge(self)
        
    def addNetwork(self,network):
        self._network = network
        
    def setExponent(self,*args):
        pass
    
    def updatefrictionFactor(self):
        pass   
 
class Status(Enum):
    pass
  
class Pipe(Edge):
    
    def __init__(self,nodeFrom,nodeTo,roughness,length,diameter):
        super().__init__(nodeFrom, nodeTo)
        ## pipe related.
        
        self.roughness = roughness 
        self.length = length 
        self.hydraulicDiameter = diameter
        self.area = None    
        
        self.reynolds = None
        self.frictionFactor = 1
        self.resistance = None
        self.gravity = 9.81
        self.k_value = 0
        self.minorResistance = 0
        
        self.calculateArea()
        self.setExponent()
        
    def setDiameter(self,diameter):
        self.hydraulicDiameter = diameter
        self.calculateArea()
         
    def calculateArea(self):
        self.area = np.pi*self.hydraulicDiameter/4
        
    def setKValue(self,k_value):
        self.k_value = k_value
    
    def setExponent(self):
        self.N_value = 2
        
    def calculateMinorLosses(self):
        self.minorResistance = self.k_value * 8 / (np.pi * self.gravity * self.hydraulicDiameter**4)   
        
    def updateReynolds(self):
        #self.reynolds = self.massFlowrate * self.hydraulicDiameter /(self._network._fluid.viscosity * self.area )
        self.reynolds = self.flowrate * self.hydraulicDiameter * self._network._fluid.density /(self._network._fluid.viscosity * self.area)
         
    def updatefrictionFactor(self):
        self.updateReynolds()
        if self.reynolds > 2000:
            # Use the colbrook Equation
            a = self.roughness/(3.7*self.hydraulicDiameter)
            b = 2.51 / (self.roughness*np.sqrt(self.frictionFactor))
            c = -2*np.log10(a + b)
            self.frictionFactor = (1/c)**2
        else:
            # Use the laminar value
            self.frictionFactor = 64/self.reynolds
            
    def updateResistance(self):
        self.updatefrictionFactor()
        self.calculateMinorLosses()
        self.resistance = (8 * self.frictionFactor * self.length) / (np.pi**2 * self.gravity *  self.hydraulicDiameter**5) + self.minorResistance
        return self.resistance       
    
class Pump(Edge): 
    pass


class Component(Edge):
    pass

    
class Fluid:
    def __init__(self,density,viscosity):
        self.density = density
        self.viscosity = viscosity
            

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
        
        # regenerate the resistances for each eage
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
        i = 0
        for e in self._edges:
            self.N[i,i] = e.N_value
            i += 1
        
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
        
        Neye = np.eye(np.size(self.N,0))
        
        # note n = 2 at the moment. This needs to be refactored to use the n matrix and N matrix function.
        schur = np.linalg.inv(self.connectivity.transpose() @ np.linalg.inv(self.G) @ self.connectivity)
        
        print("shur:\n")
        print(schur)
        
        part = self.connectivity.transpose() @ ((1-2) * self.flowVector - np.linalg.inv(self.G) @ (self.fixed @ self.fixedHeadVector)) - 2 * self.demandVector
        
        print("part:\n")
        print(part)

        self.unkownHeadVector = schur @ part
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
        
        Neye = np.eye(np.size(self.N,0))
                 
        self.flowVector = (Neye-np.linalg.inv(self.N)) @ self.flowVector + np.linalg.inv(self.N) @ np.linalg.inv(self.G) @ (self.connectivity @ self.unkownHeadVector + self.fixed @ self.fixedHeadVector)
        print(self.flowVector)
        
        i = 0
        for e in self._edges:
            e.flowrate = self.flowVector[i]
            i += 1
        
    def itterate(self,ittLimit):
        itterations = 0
        while ( itterations<ittLimit and ((self.calculateFlowResidual() > self.flowTolerance) or (self.calculatePressureResidual() > self.pressureTolerance))):
            print(f'Itterations: {itterations}')
            self.generateGMatrix()       
            self.itterateHead()
            self.itterateFlow()
            itterations += 1


