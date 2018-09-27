# File name: sim-table5.R
# 
# A code for simulation studies in "Factor-driven Two-regime Regression"
# : This code generates Table 5 in the paper and save it as "../results/table5.txt"
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
table5 = data.frame(matrix(NA, 9,6))
colnames(table5) = c('Statistic', 'Method', 'd_g=2', 'd_g=3', 'd_g=4', 'd_g=5')  # Note that we use the parameter named d.f for d_g
table5[,1] = c(rep('Min',2), rep('Median',2), rep('Mean',2), rep('Max',2), rep('Conv. Rate',1))
table5[,2] = c(rep(c('Iter(zeta=1.0)', 'Joint'),4), 'Iter(zeta=1.0)')

set.seed(325618)

result = matrix(NA, n.rep, 4)
colnames(result) = c('time.iter10', 'time.joint', 'obj.iter10', 'obj.joint') 

for (i.d.f in c(2:5)){
  d.x = 1            
  n.obs = 200
  d.f = i.d.f
  
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
    result[i,2] = end.t.joint - start.t.joint
    result[i,4] = ex_out$objval
    print(ex_out)
    
    # Iterative Method
    # Construct the grid points by equal width
    grid.gen.option = 'fixed'
    zeta.set=c(1)
    j = 1
    for (zeta in zeta.set) {
      start.t.iter = proc.time()[3]
      fixed.width = rep(zeta,d.f)
      
      grid = gen_grid(option.grid = grid.gen.option, width=fixed.width, n.total=NULL, L.grid=L.gm, U.grid=U.gm)
      
      ex_out_iter=fadtwo(y=y,x=x,f=f, method='iter', L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2, grid=grid, max.iter=max.iter)
      end.t.iter = proc.time()[3]
      result[i,j] = end.t.iter - start.t.iter
      result[i,(j+2)] = ex_out_iter$objval
      print(ex_out_iter)
      j = j + 1
    }
    
  }
  
  # save the summary statistic for iter.(zeta=1.0) for the sample size T=i.d.f
  table5[c(1,3,5,7),paste('d_g=',i.d.f,sep='')] = round(summary(result[,1])[c(1,3,4,6)],2)
  # save the summary statistic for joint for the sample size T=i.d.f
  table5[c(2,4,6,8),paste('d_g=',i.d.f,sep='')] = round(summary(result[,2])[c(1,3,4,6)],2)
  # save the convergence rate for iter.(zeta=1.0)  
  table5[9,paste('d_g=',i.d.f,sep='')] = round(mean(abs(result[,3]-result[,4])<=eta),2)
  
}

sink(file='../results/table5.txt')
print(table5)
sink()

save.image(file=paste("../results/o-sim-table5.RData",sep=''))

