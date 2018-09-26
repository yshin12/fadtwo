# File name: sim-table4.R
# 
# A code for simulation studies in "Factor-driven Two-regime Regression"
# : This code generates Table 4 in the paper and save it as "../results/table4.txt"
#
# Last Update
#   2018-09-26
#

rm(list=ls())

source('lib_fadtwo.R')

dgp1 <- function(n.obs,d.x,d.f){
  bt = rep(1,d.x)
  dt = rep(1,d.x)
  gm = c(1,rep(2/3,d.f-1))
  
  # Data dictionary: denerate x,f,eps, and y
  x = cbind(rep(1,n.obs),matrix( rnorm( (d.x-1)*n.obs, mean = 0, sd = 1 ), nrow=n.obs, ncol=(d.x-1) ))
  f = cbind(matrix(rnorm( (d.f-1)*n.obs, mean = 0, sd = 1 ), nrow=n.obs, ncol=d.f-1), rep(-1,n.obs))
  eps = rnorm( n.obs, mean = 0, sd = 0.5 )
  y = x %*% bt + x %*% dt * (f %*% gm > 0 ) + eps

  return(list(y=y, x=x, f=f, n.obs=n.obs, bt=bt, dt=dt, gm=gm))  
}

n.rep = 100
table4 = data.frame(matrix(NA, 14,6))
colnames(table4) = c('Statistic', 'Method', 'd_x=1', 'd_x=2', 'd_x=3', 'd_x=4')
table4[,1] = c(rep('Min',3), rep('Median',3), rep('Mean',3), rep('Max',3), rep('Conv. Rate',2))
table4[,2] = c(rep(c('Iter(zeta=1.0)', 'Iter(zeta=0.1)', 'Joint'),4), 'Iter(zeta=1.0)', 'Iter(zeta=0.1)')

set.seed(325618)

result = matrix(NA, n.rep, 6)
colnames(result) = c('time.iter01', 'time.iter10', 'time.joint', 'obj.iter01', 'obj.iter10', 'obj.joint') 

for (i.d.x in c(1:4)){
  d.x = i.d.x            
  n.obs = 200
  d.f = 2
  
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
  
  max.iter=1
  
  # Gurobi options: 
  params <- list(OutputFlag=1, FeasibilityTol=1e-9, MIPGap=1e-4, TimeLimit=Inf)   
  # OutputFlag: print out outputs
  # FeasibilityTol: error allowance for the inequality constraints
  # MIPGap: Stop if the gap is smaller than this bound. Default is 1e-4
  # TimeLimit: Stop if the computation time reaches this bound. Default is infinity

  for (i in (1:n.rep)){
    # DGP  
    my.dat = dgp1(n.obs, d.x, d.f)
    y=my.dat$y
    x=my.dat$x
    f=my.dat$f
    
    # Joint Method
    start.t.joint = proc.time()[3]
    ex_out=fadtwo(y=y,x=x,f=f, method='joint', L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2)
    end.t.joint = proc.time()[3]
    result[i,3] = end.t.joint - start.t.joint
    result[i,6] = ex_out$objval
    print(ex_out)
    
    # Iterative Method
    # Construct the grid points by equal width
    grid.gen.option = 'fixed'
    zeta.set=c(0.1, 1)
    j = 1
    for (zeta in zeta.set) {
      start.t.iter = proc.time()[3]
      fixed.width = rep(zeta,d.f)
      
      grid = gen_grid(option.grid = grid.gen.option, width=fixed.width, n.total=NULL, L.grid=L.gm, U.grid=U.gm)
      
      ex_out_iter=fadtwo(y=y,x=x,f=f, method='iter', L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2, grid=grid, max.iter=max.iter)
      end.t.iter = proc.time()[3]
      result[i,j] = end.t.iter - start.t.iter
      result[i,(j+3)] = ex_out_iter$objval
      print(ex_out_iter)
      j = j + 1
    }
    
  }
  
  # save the summary statistic for iter.(zeta=1.0) for the sample size T=i.d.x
  table4[c(1,4,7,10),paste('d_x=',i.d.x,sep='')] = round(summary(result[,2])[c(1,3,4,6)],2)
  # save the summary statistic for iter.(zeta=0.1) for the sample size T=i.d.x
  table4[c(2,5,8,11),paste('d_x=',i.d.x,sep='')] = round(summary(result[,1])[c(1,3,4,6)],2)
  # save the summary statistic for joint for the sample size T=i.d.x
  table4[c(3,6,9,12),paste('d_x=',i.d.x,sep='')] = round(summary(result[,3])[c(1,3,4,6)],2)
  # save the convergence rate for iter.(zeta=1.0)  
  table4[13,paste('d_x=',i.d.x,sep='')] = round(mean(abs(result[,5]-result[,6])<=eta),2)
  # save the convergence rate for iter.(zeta=0.1)  
  table4[14,paste('d_x=',i.d.x,sep='')] = round(mean(abs(result[,4]-result[,6])<=eta),2)
}

sink(file='../results/table4.txt')
print(table4)
sink()

save.image(file=paste("o-sim-table4.Rdata",sep=''))

