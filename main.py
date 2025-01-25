
from Models import Node, Pipe, Fluid, Network

# Set up the network
net0 = Network()     

fl0 = Fluid(1000,1*10**-3) 

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()
n4 = Node()

e01 = Pipe(n0,n1,1*10**-5,100,0.05)
e12 = Pipe(n1,n2,1*10**-5,100,0.05)
e13 = Pipe(n1,n3,1*10**-5,100,0.05)
e14 = Pipe(n1,n4,1*10**-5,100,0.05)

n2.setDemand(10/3600)
n3.setDemand(5/3600)
n4.setFixedPressure(10)
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