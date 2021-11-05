# Runner 2D

Optimize the running speed of your musculoskeletal model



## Installation instruction

Requirements:

- Matlab (tested in 2020b)
- Casadi API ()
- Opensim API ()



## GUI

![](App/figs/GUI_PrintScreen.png)

Overview of the graphical user interface. The goal is maximize the running speed. You have a bonus strength of 1000N to adjust the maximal isometric force of three muscles in the model (Soleus, Gastrocnemius and Rectus femoris). You can also adapt the normalized stiffness of the Achilles tendon (that is in series with the contractile elements of the Soleus and Gastrocnemius). 



## Background information

The running motion is predicted by solving an optimal control problem. We optimize the muscle excitations to maximize the running speed. We used a direct collocation approach to formulate the problem in Casadi and solve the resulting NLP in IPOPT.

Some additional information:

- Trapezoidal integration scheme
- Small term in objective function related to minimization of muscle activations squared and the controls in our simulation (to improve convergence)
- Assumed left-right symmetry which enables us to simulate half a gait cycle
- 1 contact sphere on calcaneus (heel) and two contact spheres on toes
- Used polynomials to approximate muscle-tendon-length and moment arms (see FitPolynomials)