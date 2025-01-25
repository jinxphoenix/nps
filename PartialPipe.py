from abc import ABC, abstractmethod

class PartialPipe(ABC):
    @abstractmethod
    def updateResistance(self,flowrate): 
        pass