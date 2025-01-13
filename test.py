import jax
import jax.numpy as jnp
# from jax import grad, jit, vmap
from jax import random
# from jax import jacfwd, jacrev
import matplotlib.pyplot as plt
key = random.key(0)


P1 = 5
P2 = 4
P0j = [P1]

R01 = .1
R02 = .1
R0j = [R01]


M0 = -1

# def massFunction (Pi,nodes, edges, Mi):
    
#     val = 0
    
#     for i in range(0,len(nodes),1):
    
#         val = val + ((Pi - nodes[i])/edges[i])**(1/2)
        
#     return val - M0
        

