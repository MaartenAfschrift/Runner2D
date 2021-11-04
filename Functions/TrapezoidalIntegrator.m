function error = TrapezoidalIntegrator(x,x1,xd,xd1,dt)

error = (x1-x) - (0.5*dt*(xd+xd1));