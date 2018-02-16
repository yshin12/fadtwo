# example_fadtwo.R
# 
# A simple exmple of the factor-driven two-regime regression estimation procedure.
#
# Last Update
#   2018-01-29  Add the 'iterative method' part
#   2018-01-27  Original code

rm(list=ls())
time.start = proc.time()[3]

options(digits=7)
source('lib_fadtwo_selection.R')
#sink('example_fadtwo.out')
dgp1 <- function(n.obs,d.x,d.f){
  bt = rep(1,d.x)
  dt = rep(1,d.x)
  gm = c(1,0.5417,0,0,0.5417)
  
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
d.f = 5       # 3 non-zoro coefficients (f_1, f_2, f_{d.f})
p = 3

L.bt = rep(-3,d.x)
U.bt = rep(3,d.x)

L.dt = rep(-3,d.x)
U.dt = rep(3,d.x)

L.gm1 = c(1)
U.gm1 = c(1)

L.gm2 = c(rep(-3,p))
U.gm2 = c(rep(3,p))

L.gm3 = -3
U.gm3 = 3

L.p = 0
U.p = 1

tau1 = 0.05
tau2 = 0.95

#ld = 0
ld = 2*0.5^2/n.obs               #AIC
#ld = log(n.obs)*0.5^2/n.obs     #BIC

re = 100
bt.all = rep(NA,re)
dt.all = rep(NA,re)
gm.all = matrix(NA,re,d.f)
correct.select.all = rep(NA,re)

# eta: the size of effective zero
eta = 1e-6
# FeasibilityTol: error allowance for the inequality constraints
params <- list(OutputFlag=1, FeasibilityTol=1e-9)      # Parameters for Gurobi


for (c.rep in (1:re)){
  my.dat = dgp1(n.obs, d.x, d.f)
  y=my.dat$y
  x=my.dat$x
  f=my.dat$f
  
  bt.0 = my.dat$bt           # Beta: true value from DGP1
  dt.0 = my.dat$dt           # Delta: ture value from DGP1
  gm.0 = my.dat$gm           # Gamma: true value from DGP1
  
  
  f1=f[,1]
  f2=cbind(f[,(2:(2+p-1))])
  
  #ex_out=fadtwo(y=y,x=x,f1=f1,f2=cbind(f2, rnorm(100), rnorm(100)),method='joint',L.bt=L.bt,U.bt=U.bt,L.dt=L.dt,U.dt=U.dt,L.gm1=L.gm1,U.gm1=U.gm1,L.gm2=L.gm2,U.gm2=U.gm2,L.gm3=L.gm3, U.gm3=U.gm3,
  #              L.p=L.p, U.p=U.p, tau1=tau1,tau2=tau2, eta=eta, params=params,ld=ld)
  ex_out=fadtwo(y=y,x=x,f1=f1,f2=cbind(f2),method='joint',L.bt=L.bt,U.bt=U.bt,L.dt=L.dt,U.dt=U.dt,L.gm1=L.gm1,U.gm1=U.gm1,L.gm2=L.gm2,U.gm2=U.gm2,L.gm3=L.gm3, U.gm3=U.gm3,
                L.p=L.p, U.p=U.p, tau1=tau1,tau2=tau2, eta=eta, params=params,ld=ld)
  
  print(ex_out)
  
  bt.all[c.rep] = ex_out$bt.hat
  dt.all[c.rep] = ex_out$dt.hat
  gm.all[c.rep,] = (ex_out$gm.hat) * (ex_out$gm.hat > eta )
  correct.select.all[c.rep] =  identical( (gm.all[c.rep,] != 0 ), (gm.0 !=0) )
}

save.image(file=paste('sim-sel-d.f-',d.f,'-p-',p,'.RData',sep=''))

time.end=proc.time()[3]
cpu.time = time.end-time.start

prop.correct.select = mean(correct.select.all)
cat('--------------------------------------------', '\n')
cat('Prop. correct model selection = ', prop.correct.select,'\n')
cat('--------------------------------------------', '\n')
cat('average time = ', cpu.time, '\n')
cat('--------------------------------------------', '\n')
