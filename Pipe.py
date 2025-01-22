
from PartialPipe import PartialPipe
import numpy as np

class Pipe(PartialPipe):
    
    def __init__(self,roughness,length, diameter):
        super().__init__()
        ## pipe related.
        
        self.roughness = roughness 
        self.length = length 
        self.hydraulicDiameter = diameter
        self.area = None    
        
        self.reynolds = None
        self.frictionFactor = 1
        self.resistance = None
        self.gravity = 9.81
    
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
