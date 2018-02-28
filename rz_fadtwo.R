# example_fadtwo.R
# 
# A simple exmple of the factor-driven two-regime regression estimation procedure.
#
# Last Update
#   2018-01-29  Add the 'iterative method' part
#   2018-01-27  Original code
#   2018-02-08  custimaization for the empirical model

rm(list=ls())

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

# Maximum iterations if method is 'iter'
K.bar = 2

# eta: the size of effective zero
eta = 1e-6

# Gurobi options: 
params <- list(OutputFlag=1, FeasibilityTol=1e-9, MIPGap=1e-4, TimeLimit=Inf)   
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
f = cbind(my.ts[(n.lags:(n.obs.raw-1)),c('unemp')],-1)
#f = cbind(my.ts[(n.lags:(n.obs.raw-1)),c('unemp', 'tbill')],-1)
#f = cbind(my.ts[(n.lags:(n.obs.raw-1)),'unemp'], my.ts[(n.lags:(n.obs.raw-1)),'tbill'], my.ts[(n.lags:(n.obs.raw-1)),'y'],-1)

# Number of factors, dim(f)
d.f = ncol(f)

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

# Joint Estimation Algorithm
if (method == 'joint') {
  time.start = proc.time()[3]

  est.out=fadtwo(y=y,x=x,f=f, method='joint', L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2)
  print(est.out)
  
  time.end = proc.time()[3]
  runtime = time.end - time.start
  cat('Runtime     =', runtime, 'sec','\n')  
}

# Iterative Estimation Algorithm
if (method == 'iter'){
  time.start = proc.time()[3]
  
  # Generate the grid points
  zeta=2
  grid.base = seq(from=-Bnd.Const, to=Bnd.Const, by =zeta)
  if (d.f == 2){
    grid=expand.grid(1,grid.base)  
  } else if (d.f == 3){
    grid=expand.grid(1,grid.base, grid.base)  
  } else if (d.f == 4 ){
    grid=expand.grid(1,grid.base, grid.base, grid.base)  
  } else if (d.f == 5) {
    grid=expand.grid(1,grid.base, grid.base, grid.base, grid.base)  
  } else {
    # Define your grid based on the number of factors in your model. For example,
    # if d.f=6, then
    # grid=expand.grid(1,grid.base, grid.base, grid.base, grid.base, grid.base)  
  }
  
  
  est.out=fadtwo(y=y,x=x,f=f, method='iter', L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2, grid=grid, max.iter=K.bar)
  print(est.out)

  time.end = proc.time()[3]
  runtime = time.end - time.start
  cat('Runtime     =', runtime, 'sec','\n')  
}


#-----------------------------------------------------------------------------------------------------------
#
# Post-estimation analysis
#
#----------------------------------------------------------------------------------------------------------
# Save the esitmates
bt.hat = est.out$bt.hat
dt.hat = est.out$dt.hat
gm.hat = est.out$gm.hat
objval = est.out$objval

# Estimate the impulse-reponse function using the local projection method
# y (output) equation
state.y = (f %*% gm.hat >= eta)
x.nl.y = cbind(x*as.numeric(1-state.y), x*as.numeric(state.y))
colnames(x.nl.y)=c(paste(colnames(x),'.state1',sep=''),paste(colnames(x),'state2',sep='.'))  # state 2: f'gm > 0 

# Columns of shocks are subtracted
x.subt.shock.y = x.nl.y[,!(colnames(x.nl.y)=='.newsy.state1' | (colnames(x.nl.y)=='.newsy.state2'))]
x.shock.y = x.nl.y[,c('.newsy.state1','.newsy.state2')]
IRF.nl.y=irf_local_projection(y=y, x=x.subt.shock.y, s=x.shock.y, h=20, print=FALSE)$irf

# g (gov't expenditure) equation
m.g = lm(g~x.nl.y-1)
bt.hat.g = coef(m.g)[1:d.x]
dt.hat.g = coef(m.g)[-c(1:d.x)]
IRF.nl.g = irf_local_projection(y=g, x=x.subt.shock.y, s=x.shock.y, h=20, print=FALSE)$irf

# Draw IRFs
pdf(file="irf-y.pdf")
par(cex.lab=1.5, mar=c(5,6,4,1)+.1) 
plot(c(0:20),IRF.nl.y[,1], type='l', col='blue', lty=2, ylab='GDP', xlab='Quarter', xlim=c(0,20), xaxp=c(0,20,4),  ylim=c(min(c(IRF.nl.g[,1],IRF.nl.g[,2])),max(c(IRF.nl.g[,1],IRF.nl.g[,2]))))
abline(h=0,lty=1) 
lines(c(0:20),IRF.nl.y[,2], type='l', col='red')
legend(0, 1.0, c('S-1 (f*gm =< 0)', 'S-2 (f*gm > 0)'), lty=c(2,1), col = c('blue','red'))
dev.off()

pdf(file="irf-g.pdf")
par(cex.lab=1.5, mar=c(5,6,4,1)+.1) 
plot(c(0:20),IRF.nl.g[,1], type='l', col='blue', lty=2, ylab='Government Spending', xlab='Quarter', xlim=c(0,20), xaxp=c(0,20,4), ylim=c(min(c(IRF.nl.g[,1],IRF.nl.g[,2])),max(c(IRF.nl.g[,1],IRF.nl.g[,2]))))
abline(h=0,lty=1) 
lines(c(0:20),IRF.nl.g[,2], type='l', col='red')
legend(0, 1.0, c('S-1 (f*gm =< 0)', 'S-2 (f*gm > 0)'), lty=c(2,1), col = c('blue','red'))
dev.off()

# Draw recession periods on the output series
pdf(file="recessions.pdf")
recess = ggplot(data.frame(my.ts, dat.year))+
  geom_line(mapping=aes_string(x="dat.year", y="y"))+
  annotate("rect", fill = "blue", alpha = 0.2, 
           xmin = dat.year[state.y], xmax = dat.year[state.y]+0.25,
           ymin = -Inf, ymax = Inf) +
  xlab('Year') +
  ylab('GDP')
print(recess)
dev.off()


# Generate new data set for STATA - IV reg
dat.temp=read.csv(file='rz_dat_original.csv', header=T)
fstate = rep(NA,nrow(dat.temp))
n.miss = length(fstate)-length(state.y)
fstate[-(1:n.miss)] = as.integer(state.y)
dat.temp = cbind(dat.temp,fstate)
write.csv(dat.temp, file='rz_dat_updated.csv', na='')




