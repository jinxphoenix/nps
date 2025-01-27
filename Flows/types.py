
from enum import Enum


class Diameter:
    def __init__(self,diameter):
        self.diameter = diameter
        
    def __getattribute__(self, name):
        return self.__dict__[f"_{name}"]
    
    def __setattr__(self, name, value):
        self.__dict__[f"_{name}"] = value
        
class PipeStatus(Enum):
    ENABLED = 1
    DISABLED = 2
    
    