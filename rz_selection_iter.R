# example_fadtwo.R
# 
# A simple exmple of the factor-driven two-regime regression estimation procedure.
#
# Last Update
#   2018-01-29  Add the 'iterative method' part
#   2018-01-27  Original code

rm(list=ls())

options(digits=7)

library('ggplot2')
source('lib_fadtwo.R')


# ------------------------------------------------------------------------------------
# 
# Set the environment variables
#
#-------------------------------------------------------------------------------------
# Choose the data length
# 'short': 1920:1 - 2015:4
# 'long' : 1889:1 - 2015:4
# Notice that the 'tbill' series starts from 1920:1 and you should choose 'short' for
# using tbill as a factor
dat.length = 'short' # 'short' or 'long'

# Estimation algorithm
method = 'iter' # 'joint' or 'iter'

# eta: the size of effective zero
eta = 1e-6

# Gurobi options: 
params <- list(OutputFlag=1, FeasibilityTol=1e-9, MIPGap=10e-4, TimeLimit=Inf)   
# OutputFlag: print out outputs
# FeasibilityTol: error allowance for the inequality constraints
# MIPGap: Stop if the gap is smaller than this bound. Default is 1e-4
# TimeLimit: Stop if the computation time reaches this bound. Default is infinity


# Minimum and maximum bounds for the propotion of each regime
tau1 = 0.05
tau2 = 0.95

# Number of lags in the regressor part, x
n.lags=4                                    

# Number of time horizon for the impulse-response function
h=20                                    
# ------------------------------------------------------------------------------------
# 
# End of seeting the environment variables
#
#-------------------------------------------------------------------------------------


# Read Data
raw.ts=read.csv(file='dat_rz.csv', header=TRUE)
if (dat.length == 'short'){
  my.ts=ts(raw.ts[-c(1:124),], start=c(1920,1), end=c(2015,4), frequency=4)
  dat.year =  time(my.ts)
  my.ts=my.ts[,c('y', 'g', 'newsy', 'unemp', 'tbill')]
}
if (dat.length == 'long'){
  my.ts=ts(raw.ts, start=c(1889,1), end=c(2015,4), frequency=4) 
  dat.year =  time(my.ts)[-(1:4)]
  my.ts=my.ts[-(1:4),]  # Missing observations in the first 4 quarters. Dropped.  
  my.ts=my.ts[,c('y', 'g', 'newsy', 'unemp')]
}
n.obs.raw=nrow(my.ts)

# output (gap)
y = my.ts[-(1:n.lags),'y']     
# govn't expenditure
g = my.ts[-(1:n.lags),'g']     
x = cbind(1, my.ts[-(1:n.lags),3], my.ts[(n.lags:(n.obs.raw-1)),c(1:3)], my.ts[((n.lags-1):(n.obs.raw-2)),c(1:3)], my.ts[((n.lags-2):(n.obs.raw-3)),c(1:3)], 
          my.ts[((n.lags-3):(n.obs.raw-4)),c(1:3)] )
colnames(x) = c('.const','.newsy','.L1.y','.L1.g','.L1.newsy','.L2.y','.L2.g','.L2.newsy','.L3.y','.L3.g','.L3.newsy','.L4.y','.L4.g','.L4.newsy')
# Rearrange the regressors so that they are ordered by lags of each regressors
x = x[,c('.const','.newsy',paste('.L',(1:4),'.newsy',sep=''),paste('.L',(1:4),'.y',sep=''),paste('.L',(1:4),'.g',sep=''))]   
# Number of observations
n.obs=nrow(x)
# Number of covariates, dim(x)
d.x = ncol(x)
#-----------------------------------------------------------------------------------------------------------------------------
#
# Choose factors. Again, notice that dat.length should be 'short' to add 'tbill' as a facotr
# When you change the factor design, you also need to change grid, so that the dimension of grid matches the number of factors
#
#-----------------------------------------------------------------------------------------------------------------------------
f = cbind(my.ts[(n.lags:(n.obs.raw-1)),c('unemp','tbill')],-1)
#f = cbind(my.ts[(n.lags:(n.obs.raw-1)),'unemp'], my.ts[(n.lags:(n.obs.raw-1)),'tbill'], my.ts[(n.lags:(n.obs.raw-1)),'y'],-1)

# Number of factors, dim(f)
d.f = ncol(f)
f1=as.matrix(f[,1])
f2=as.matrix(f[,2])
p = ncol(f2)   # Selectino of T-bill given unemp is included in the model

#-----------------------------------------------------------------------------------------------------------------------------
#
# Choose the size of the parameter sets
#
#-----------------------------------------------------------------------------------------------------------------------------
Bnd.Const = 20
L.bt = rep(-Bnd.Const,d.x)
U.bt = rep(Bnd.Const,d.x)

L.dt = rep(-Bnd.Const,d.x)
U.dt = rep(Bnd.Const,d.x)

L.gm = c(1,rep(-Bnd.Const,d.f-1))
U.gm = c(1,rep(Bnd.Const,d.f-1))



L.gm1 = c(1)
U.gm1 = c(1)

L.gm2 = c(rep(-Bnd.Const,p))
U.gm2 = c(rep(Bnd.Const,p))

L.gm3 = -Bnd.Const
U.gm3 = Bnd.Const

zeta=2
grid.base = seq(from=L.gm2, to=U.gm2, by =zeta)
grid=expand.grid(1,grid.base,grid.base)

L.p = 0
U.p = 1

ld = log(n.obs)*var(y)/n.obs     #BIC

# K.bar: number of maximum iterations
K.bar = 1


# eta: the size of effective zero
eta = 1e-6

ex_out=fadtwo_selection(y=y, x=x, f1=f1, f2=f2, method=method, L.gm1=L.gm1, U.gm1=U.gm1, L.gm2=L.gm2, U.gm2=U.gm2, L.gm3=L.gm3, U.gm3=U.gm3, L.p=L.p, U.p=U.p, tau1=tau1, tau2=tau2, eta=eta, params=params, grid=grid, max.iter=K.bar, ld=ld, p=p)

print(ex_out)

# Iterative Method

#sink()
