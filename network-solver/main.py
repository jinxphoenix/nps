
from Models import Node, Pipe, Fluid, Network, Pump

# Set up the network
net0 = Network()     

fl0 = Fluid(1000,1*10**-3) 

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()
# n4 = Node()
# n5 = Node()

e01 = Pipe(n0,n1,1*10**-5,100,0.1)
# e12 = Pipe(n1,n2,1*10**-5,100,0.05)
e12 = Pump(n1,n2,-6680.879290080472,49.93200286913462,1.493199982498533)
e23 = Pipe(n2,n3,1*10**-5,100,0.1)
# e14 = Pipe(n1,n4,1*10**-5,100,0.05)
# e45 = Pipe(n4,n5,1*10**-5,100,0.05)

# n2.setDemand(20/3600)
# n3.setDemand(80/3600)
n0.setFixedPressure(20)
n3.setFixedPressure(30)
# n5.setFixedPressure(5)
# e14.set_k_value(2)
# n0.setHeight(5)
n3.setHeight(20)

net0.addFluid(fl0)
net0.addNode(n0)
net0.addNode(n1)
net0.addNode(n2)
net0.addNode(n3)
# net0.addNode(n4)
# net0.addNode(n5)

net0.addEdge(e01)
net0.addEdge(e12)
net0.addEdge(e23)
# net0.addEdge(e13)
# net0.addEdge(e14)
# net0.addEdge(e45)

# initialise and run the simulation

net0.initialise()
net0.itterate(10)