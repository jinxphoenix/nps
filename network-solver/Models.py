
from enum import Enum
import numpy as np
# from Flows.types import PipeStatus

class Fluid:
    def __init__(self,density,viscosity):
        self.density = density
        self.viscosity = viscosity

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
        
        # Basic edge parameters
        self.ID = Edge._nextID
        self._nodeFrom: Node = nodeFrom
        self._nodeTo: Node = nodeTo
        Edge._nextID += 1 
        self._network: Network = None
        
        # utility 
        self.status = None

        # basic expected properties
        self.flowrate = 0.1 # dmmy vale for development
        self.previousFlowrate = self.flowrate - 1 # Initiallised to a value of not the prev flow rate so first itteration works
        self.exponent = None
        self.resistance_offset = None
        self.resistance = None  
                  
    def attachFromNode(self):
        self._nodeFrom.attachEdge(self)
        
    def attachToNode(self):
        self._nodeTo.attachEdge(self)
        
    def addNetwork(self,network):
        self._network = network
        
    def set_exponent(self,N):
        pass
    
    def get_resistance_offset(self):
        pass
    
    def updateResistance(self):
        pass   
 
class Status(Enum):
    pass

class PartialPipe(Edge):
    def __init__(self, nodeFrom, nodeTo):
        super().__init__(nodeFrom, nodeTo)
        
        # nominal
        self.diameter = None
        self.area = None
        self.velocity = None
        self.reynolds = None
        
        # inlet
        self.inletDiameter = None
        self.inletArea = None
        self.inletVelocity = None
        self.inletReynolds = None
        
        # outlet
        self.outletDiameter = None
        self.outletArea = None
        self.outletVelocity = None
        self.outletReynolds = None
        
        # utilities values
        self.gravity = 9.81
        
    @staticmethod    
    def calculateArea(diameter: float):
        return np.pi /4 * diameter**2 
    
    @staticmethod
    def calculateVelocity(flowrate,area):
        return flowrate / area
    
    @staticmethod
    def calculateReynolds(fluid: Fluid, flowrate, charLength, area):
        """Calculates the reynolds number for a full pipe.

        Args:
            fluid (Fluid): Object from network containing fluid properties
            flowrate (_type_): flow rate in m3/s
            charLength (_type_): characteristic length (i.e diameter)
            area (_type_): cross sectional area

        Returns:
            _type_: reynolds number (unitless)
        """
        return flowrate * charLength * fluid.density /(fluid.viscosity * area)

    def set_exponent(self, N):
        self.exponent = N
    
class Pipe(PartialPipe):
    
    def __init__(self,nodeFrom,nodeTo,roughness,length,diameter):
        super().__init__(nodeFrom, nodeTo)
        ## pipe related.
        
        self.set_diameter(diameter)
        self.length = length   
        self.roughness = roughness 
        self.k_value = 0
        
        # emergent from flow conditions
        self.frictionFactor = 1     
        self.minorResistance = 0
        
        # set teh exponent to two as hard coded using darcey for pipes
        self.set_exponent(2)
        
        # this is a pipe, so it has no resistance offset at zero flow
        self.resistance_offset = 0
        
    def set_diameter(self,diameter):
        self.diameter = diameter
        self.inletDiameter = diameter
        self.outletDiameter = diameter
        
        # update areas when the diameter is changed.
        self.set_area()
         
    def set_area(self):
        self.area = self.calculateArea(self.diameter)
        self.inletArea = self.calculateArea(self.inletDiameter)
        self.outletArea = self.calculateArea(self.outletDiameter)
        
    def set_k_value(self,k_value):
        self.k_value = k_value
        
    def calculateMinorLosses(self):
        self.minorResistance = self.k_value * 8 / (self.gravity * np.pi**2 * self.diameter**4)   
        
    def update_velocities(self):
        pass
    
    def updateReynolds(self):      
        self.reynolds = self.calculateReynolds(self._network._fluid,self.flowrate,self.diameter,self.area)
        self.inletReynolds = self.calculateReynolds(self._network._fluid,self.flowrate,self.inletDiameter,self.inletArea)
        self.outletReynolds = self.calculateReynolds(self._network._fluid,self.flowrate,self.outletDiameter,self.outletDiameter)
         
    def updatefrictionFactor(self):
        self.updateReynolds()
        if self.reynolds > 2000:
            # Use the colbrook Equation
            a = self.roughness/(3.7*self.diameter)
            b = 2.51 / (self.roughness*np.sqrt(self.frictionFactor))
            c = -2*np.log10(a + b)
            self.frictionFactor = (1/c)**2
        else:
            # Use the laminar value
            self.frictionFactor = 64/self.reynolds
    
    def get_resistance_offset(self):
        return self.resistance_offset
            
    def updateResistance(self):
        self.updatefrictionFactor()
        self.calculateMinorLosses()
        self.resistance = (8 * self.frictionFactor * self.length) / (np.pi**2 * self.gravity *  self.diameter**5) + self.minorResistance
        return self.resistance       
   
class Pump(PartialPipe):
    def __init__(self, nodeFrom, nodeTo, alpha, beta, n):
        super().__init__(nodeFrom, nodeTo)

        # pump related
        self.set_exponent(n)
        self.set_flow_coefficient(alpha)
        self.set_resistance_offset(beta)
 
    def set_inlet_diameter(self,diameter):
        self.inletDiameter = diameter
        self.inletArea = self.calculateArea(diameter)
        
    def set_outlet_diameter(self,diameter):
        self.outletDiameter = diameter
        self.outletArea = self.calculateArea(diameter)
    
    def set_resistance_offset(self,value):
        self.resistance_offset = value
    
    def get_resistance_offset(self):
        return -self.resistance_offset
    
    def set_flow_coefficient(self,value):
        self.flow_coefficient = value
    
    def updateResistance(self):
        return -self.flow_coefficient
        
class Component(Edge):
    pass
            
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
                
    def addEdge(self,edge: Edge):
        if edge not in self._edges:
            self._edges.append(edge)
            edge.addNetwork(self)
            
    def addNode(self,node: Node):
        if node not in self._nodes:
            self._nodes.append(node)
            if node.fixedPressure == False:
                self._unodes.append(node)
            else:
                self._fnodes.append(node)
    
    def addFluid(self,fluid):
        self._fluid = fluid
                    
    def generateConnectivityMatrix(self):
        # this is the A12 matrix num edges x bum unkown pressure nodes
        self.connectivity = np.zeros((len(self._unodes),len(self._edges)))
        i = 0
        for n in self._unodes:
            j = 0
            for e in self._edges:
                if e._nodeFrom == n:
                    self.connectivity[i, j] = -1
                                        
                elif e._nodeTo == n:
                    self.connectivity[i, j] = 1
                j += 1
            i += 1
        self.connectivity = np.transpose(self.connectivity)    
                
    def generateFixedMatrix(self):
        # this is the A10 matri
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
        self.fixed = np.transpose(self.fixed)
      
    def generateResistanceMatrix(self):
        self.resistances = np.zeros((len(self._edges),len(self._edges)))
        i = 0
        for e in self._edges:
            self.resistances[i,i] = e.updateResistance()
            i += 1
    
    def generateResistanceOffsetMatrix(self):
        # this creates a matrix of beta's, divided by the right amount of Q
        pass
        self.resistance_offsets = np.zeros((len(self._edges),len(self._edges)))
        i = 0
        for e in self._edges:
            self.resistance_offsets[i,i] = e.get_resistance_offset()
            i += 1
    
    def createUnkownHeadVector(self):
        """Create a vector of unkown head nodes
        for each node in the unkown heads list, get the nodes pressure, and add the height to the node.
        """
        self.unkownHeadVector = np.zeros((len(self._unodes),1))
        i = 0
        for n in self._unodes:
            self.unkownHeadVector[i] = n.pressure + n.height
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
    
    def generateA11Matrix(self):
        # this is the A11 matrix, which is the resistances multipled by the prev itteration flows
        self.generateResistanceMatrix()
        self.generateResistanceOffsetMatrix()
        
        Qabs = np.abs(self.flowVector)
        Npow = np.diag(self.N)[:,np.newaxis]
        Qpow = np.power(Qabs,(Npow-1))
        A11naked = np.multiply(self.resistances, Qpow[:])
        B11 = np.multiply(self.resistance_offsets, (1/self.flowVector)[:])
        
        # create the diagonal matrix by multiplying row wise, the resistance matrix with the abs flow vector
        self.A11 = A11naked + B11
        # print(self.A11)
        
    def generateGMatrix(self):
        # this is the G matrix (it is A11 frotn matrix multiplied with the N matrix.)
        self.generateA11Matrix()
        
        # create the A11* matrix
        A11_star = self.A11 - np.multiply(self.resistance_offsets, (1/self.flowVector)[:])
        
        # then create the G matrix
        self.G = self.N @ A11_star
        # print(f'G vector:\n{self.G}')
        
    def generateNMatrix(self):
        """ generates an N matrix based on the exponents required in the Jacobian for each edge
            For each edges, it examines what exponent is used for flow (darcey = 2, hazen = 1.852)
            and creates a diagonal matrix for each edge with that value of N.
            
            This is used as part of the jacobian (first derivative) and is thus used by 
        """
        self.N = np.eye(len(self._edges),len(self._edges))
        i = 0
        for e in self._edges:
            self.N[i,i] = e.exponent
            i += 1
        
    def createFixedHeadVector(self):
        """Create a vector of the fixed heads in the network
        For each node in the fixed head node list, add to the fixed head node vector, and add the grade height of node to teh pressure
        """
        self.fixedHeadVector = np.zeros((len(self._fnodes),1))
        i = 0
        for n in self._fnodes:
            self.fixedHeadVector[i] = n.pressure + n.height
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
        
    def setInitalConditions(self):
        """Sets all flow and pressure values to values sutiable to start calculation
        """
        pass    
                
    def initialise(self):
        """Initialises all the matrices and vectors required for the calculation.
        """
        
        self.generateConnectivityMatrix()
        self.generateFixedMatrix()
        self.generateNMatrix()     
        self.createDemandVector()
        self.createFixedHeadVector()
        
        self.createFlowVector()
        self.createUnkownHeadVector()
        
    
    def itterateHead(self):
        # itterate the head value
        
        # first store the previous itterations pressures
        i = 0
        for n in self._unodes:
            n.previousPressure = n.pressure
            i += 1     
                    
        # note n = 2 at the moment. This needs to be refactored to use the n matrix and N matrix function.
        A = self.connectivity.transpose() @ np.linalg.inv(self.G) @ self.connectivity
        b = -( self.connectivity.transpose() @ np.linalg.inv(self.G) @ (self.A11 @ self.flowVector + self.fixed @ self.fixedHeadVector) -- (self.connectivity.transpose() @ self.flowVector - self.demandVector)  )
               
        # print(f"A:\n {A}")
        # print(f"b:\n {b}")

        self.unkownHeadVector = np.linalg.inv(A) @ b 
        print(f'unkown Head vector:\n{self.unkownHeadVector}')
        
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
        
        # create I matrix same size as N
        Neye = np.eye(np.size(self.N,0))
                 
        self.flowVector = (Neye - np.linalg.inv(self.G) @ self.A11) @ self.flowVector - np.linalg.inv(self.G) @ (self.connectivity @ self.unkownHeadVector + self.fixed @ self.fixedHeadVector)
        
        print(f'flow vector:\n{self.flowVector}')
        
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
