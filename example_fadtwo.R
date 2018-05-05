# example_fadtwo.R
# 
# A simple exmple of the factor-driven two-regime regression estimation procedure.
#
# Last Update
#   2018-01-29  Add the 'iterative method' part
#   2018-01-27  Original code

rm(list=ls())

source('lib_fadtwo.R')

dgp1 <- function(n.obs,d.x,d.f){
  bt = rep(1,d.x)
  dt = rep(1,d.x)
  gm = c(1,rep(0.5716,d.f-1))
  
  # Data dictionary: denerate x,f,eps, and y
  x = cbind(rep(1,n.obs),matrix( rnorm( (d.x-1)*n.obs, mean = 0, sd = 1 ), nrow=n.obs, ncol=(d.x-1) ))
  f = cbind(matrix(rnorm( (d.f-1)*n.obs, mean = 0, sd = 1 ), nrow=n.obs, ncol=d.f-1), rep(-1,n.obs))
  eps = rnorm( n.obs, mean = 0, sd = 0.5 )
  y = x %*% bt + x %*% dt * (f %*% gm > 0 ) + eps

  return(list(y=y, x=x, f=f, n.obs=n.obs, bt=bt, dt=dt, gm=gm))  
}

set.seed(325618)
n.obs = 100                 
d.x = 1
d.f = 3

L.bt = rep(-3,d.x)
U.bt = rep(3,d.x)

L.dt = rep(-3,d.x)
U.dt = rep(3,d.x)

L.gm = c(1,rep(-3,d.f-1))
U.gm = c(1,rep(3,d.f-1))

tau1 = 0.05
tau2 = 0.95

# eta: the size of effective zero
eta = 1e-6

# Gurobi options: 
params <- list(OutputFlag=1, FeasibilityTol=1e-9, MIPGap=1e-4, TimeLimit=Inf)   
# OutputFlag: print out outputs
# FeasibilityTol: error allowance for the inequality constraints
# MIPGap: Stop if the gap is smaller than this bound. Default is 1e-4
# TimeLimit: Stop if the computation time reaches this bound. Default is infinity

my.dat = dgp1(n.obs, d.x, d.f)
y=my.dat$y
x=my.dat$x
f=my.dat$f

ex_out=fadtwo(y=y,x=x,f=f, method='joint', L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2)
print(ex_out)

# Iterative Method

# Construct the grid points by equal width
grid.gen.option = 'fixed'
zeta = 0.1
fixed.width = rep(zeta,d.f)

grid = gen_grid(option.grid = grid.gen.option, width=fixed.width, n.total=NULL, L.grid=L.gm, U.grid=U.gm)

# Construct the grid random grid points
# grid.gen.option = 'random'
# total.grid.points = 1000
#grid = gen_grid(option.grid = grid.gen.option, width=NULL, n.total=total.grid.points, L.grid=L.gm, U.grid=U.gm)


ex_out_iter=fadtwo(y=y,x=x,f=f, method='iter', L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2, grid=grid)
print(ex_out_iter)
