## 2D predictive simulations





Repository for predictive simulations with 2D model.


### workflow

Create all required inputs for simulation based on opensim model using the scripts in the folder ConvertOsimModel. This script uses the following functions: 

- FitPolynomials: function used to compute polynomial coefficient that can describe the muscle geometry of the opensim model.
- CreateCasdiFunctions: function used to create all casadifunctions used to formulate the optimization problem.


### Models

Based on Gait18

- Gait18_Visual.osim: default gait18 model with contact spheres
- Gait18_LongerHamstring: default gait 18 with optimal length of hamstring adapted to 0.11 (from 0.10765).
- Gait18_shorterHamstrings: default gait 18 with slack length of hamstrings adapted to 0.3 (from 0.32)

### Noisy ankle

Implementation with additional torque on ankle joint. This is now modeled as a since wave (sin(ti)) with ti the discrete time.

$$ \tau_l - 0.1sin(5*ti) = \tau_{mus}$$

symmetry is solved as:

```matlab
% ensure that a*sin(b*N) = a*sin(0)
% this is the case if sin(b*N) = 0 and therefore if b*N = x 2*pi. with x an
% "geheel getal".
% b = (x * 2*pi)/N/. lets say x = 10
nPeriods = 5;
b = (nPeriods*2*pi)./N;
a = 10; % magnitude of the noise
xV = 0:1:N;
NoiseL = a*sin(b*xV);
NoiseR = a*sin(b*xV);

```

### Noisy ankle multiple phases



We simulate n-gait cycles with n-different sinus profiles as disturbance. To impose a similar control structure between the n-different gait cycles I added an additional term to the objective function that mimimizes the squared different between the average excitation (average of the n-cycles) and the excitation in the nth gait cycle for each muscle and on each descrete time point.



Run this example with DeafultExample_trpezoidal_NoisyAnkle_nphase which used the function f_PredSim_2D_trapezoidal_NoisyAnkle_nPhase.





