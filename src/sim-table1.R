# file name: sim-table1.R
# 
# A code for simulation studies in "Factor-driven Two-regime Regression"
# : This code generates Table 1 in the paper and save it as "../results/table1.txt"
#
# Last Update
#   2018-09-26
#
rm(list=ls())
time.start = proc.time()[3]
n.sim='001'
set.seed(n.sim)

source('lib_fadtwo.R')

dgp2 <- function(n.obs,n.crs, d.x, rho.e, rho.g, rho.x, K, sg_eps, bt, dt, phi){

  # Generating the Factor processes
    # (1) e_it = rho_i * e_it-1 + eps_it
    eps_e = matrix(rnorm(n.crs*(n.obs+1)), n.crs, n.obs+1)
    e = matrix(NA, n.crs, n.obs)
    e = cbind(rep(0,n.crs), e)    # Add the inital value
    for (i in c(2:(n.obs+1))){
      e[,i] = diag(rho.e) %*% e[,i-1] + eps_e[,i]
    }
    e = e[,c(2:(n.obs+1))]
    # (2) f_jt = alpha_j * e_it-1 + eps_it
    u = matrix(rnorm(K*(n.obs+1)), K, n.obs+1)
    g1 = matrix(NA, K, n.obs)
    g1 = cbind(rep(0,K), g1)    # Add the inital value
    for (i in c(2:(n.obs+1))){
        g1[,i] = diag(rho.g) %*% g1[,i-1] + u[,i]  
    }
    g1 = g1[,c(2:(n.obs+1))]
    # (3) X = Lambda * g1 + \sqrt{K}*e
    Lambda = matrix(rnorm(n.crs*K)*sqrt(K), n.crs, K)
    Xit = Lambda %*% g1 + sqrt(K) * e
  
  # Generate x from AR(1) with the coefficient rho.x
  x = rep(1,n.obs)
  if (d.x > 1) {
    for ( i in c(1:(d.x-1))){
      x = cbind(x, arima.sim(model=list(ar=rho.x[i]), n=n.obs)) 
    }
  }
  colnames(x)=paste('x',c(1:d.x),sep='')

  # Generate y with the true factors
  eps = rnorm( n.obs, mean = 0, sd = sg_eps )
  g = cbind(t(g1), rep(-1,n.obs))
  state = (g %*% phi > 0 )
  y = x %*% bt + x %*% dt * state + eps

  return(list(y=y, x=x, Xti=t(Xit), g=g, state0=state, Lambda=Lambda))  
}

# parameters for DGP
n.obs = 200     
n.crs = 200
n.rep = 1000
d.x = 2             # Dimension of x_t
rho_x = 0.5         # AR(1) coefficient for x_t
K = 3               # The number of non-constant factors
sg_eps = 0.5        # Standard deviation of epsilon (the error term of the model)

p = 2
feasible.f = 2+p  # This dimension includes the constant term (f1,f2,f3,-1)
rho_e = runif(n.crs, min=0.3, max=0.5)
rho_g = runif(K, min=0.2, max=0.8)

bt0 = rep(1,d.x)
dt0 = rep(1,d.x)
ap0 = c(bt0, dt0)
phi0 = c(1,2/3,0,2/3)


bnd = 3

L.bt = rep(-bnd,d.x)
U.bt = rep(bnd,d.x)

L.dt = rep(-bnd,d.x)
U.dt = rep(bnd,d.x)

L.gm = c(1,rep(-bnd,K-1))
U.gm = c(1,rep(bnd,K-1))

L.gm.full = c(1,rep(-bnd, K))
U.gm.full = c(1,rep(bnd, K))

tau1 = 0.05
tau2 = 0.95

L.gm1 = c(1)
U.gm1 = c(1)

L.gm2 = c(rep(-bnd,p))
U.gm2 = c(rep(bnd,p))

L.gm3 = -bnd
U.gm3 = bnd

L.p = 0
U.p = p

#BIC
ld = sg_eps^2*log(n.obs)/n.obs   # sig_eps^2 *log(n.obs) / n.obs


# eta: the size of effective zero
eta = 1e-6

# Gurobi options: 
params <- list(OutputFlag=1, FeasibilityTol=1e-9, MIPGap=1e-4, TimeLimit=Inf)   
# OutputFlag: print out outputs
# FeasibilityTol: error allowance for the inequality constraints
# MIPGap: Stop if the gap is smaller than this bound. Default is 1e-4
# TimeLimit: Stop if the computation time reaches this bound. Default is infinity

data.base.sim = array(NA, dim=c(n.obs,1+d.x+feasible.f,n.rep))
  
result.o = matrix(NA, n.rep, d.x+d.x+n.obs)         # bt.hat; dt.hat; state0
result.k1 = matrix(NA, n.rep, d.x+d.x+K+1+n.obs)     # bt.hat; dt.hat; gm.hat; objval; state.hat
result.k2 = matrix(NA, n.rep, d.x+d.x+feasible.f+1+n.obs)     # bt.hat(d.x); dt.hat(d.x); gm.hat(K); objval(1); state.hat(T)
result.u = matrix(NA, n.rep, d.x+d.x+2*feasible.f+1+n.obs)     # bt.hat(d.x); dt.hat(d.x); gm.hat(K); objval(1); state.hat(T)
colnames(result.o) = c(paste('bt',(1:d.x),sep=''), paste('dt',(1:d.x),sep=''), paste('state',(1:n.obs),sep=''))
colnames(result.k1) = c(paste('bt',(1:d.x),sep=''), paste('dt',(1:d.x),sep=''), paste('gm',(1:K),sep=''), 'objval', paste('state',(1:n.obs),sep=''))
colnames(result.k2) = c(paste('bt',(1:d.x),sep=''), paste('dt',(1:d.x),sep=''), paste('gm',(1:feasible.f),sep=''), 'objval', paste('state',(1:n.obs),sep=''))
colnames(result.u) = c(paste('bt',(1:d.x),sep=''), paste('dt',(1:d.x),sep=''), paste('gm',(1:feasible.f),sep=''), paste('gm0_',(1:feasible.f),sep=''), 'objval', paste('state',(1:n.obs),sep=''))

cover.o = matrix(NA, n.rep, d.x*2)
cover1.k1 = matrix(NA, n.rep, d.x*2)
cover2.k1 = matrix(NA, n.rep, d.x*2) 
cover1.k2 = matrix(NA, n.rep, d.x*2)
cover2.k2 = matrix(NA, n.rep, d.x*2) 
cover1.u = matrix(NA, n.rep, d.x*2)
cover2.u = matrix(NA, n.rep, d.x*2) 

i.rep=1
while (i.rep <= n.rep) {
  flag.dgp.issue = FALSE
  cat('----------------------- \n', i.rep, 'iterations','\n','-------------------------','\n')
  my.dat = dgp2(n.obs,n.crs, d.x, rho_e, rho_g, rho_x, K, sg_eps, bt=bt0, dt=dt0, phi=phi0) 
  y=my.dat$y
  x=my.dat$x
  Xti=my.dat$Xti
  g=my.dat$g
  state0=my.dat$state0
  J0 = (phi0 != 0)
  g0=(my.dat$g)[, J0]
  
  data.base.sim[ , ,i.rep] = cbind(y, x, g)

  # -----------------------------------------------------------------------------------------------------------------
  # O: Estimate the model with true states (Oracle)
  # -----------------------------------------------------------------------------------------------------------------
    cat('---------------------------------------- \n')
    cat('O: known states \n')
    cat('---------------------------------------- \n \n')
    reg = cbind(x, x*matrix(state0,n.obs,d.x))
    colnames(reg) = c(paste('x',c(1:d.x),sep=''), paste('x*state',c(1:d.x),sep=''))
    out_o1 = lm(y~reg-1)
    print(out_o1)
    result.o[i.rep,] = c(coef(out_o1), state0)
    ci.o = confint(out_o1)
    cover.o[i.rep,] = (ci.o[,1] <= c(bt0,dt0)) & (ci.o[,2] >= c(bt0,dt0))

  # -----------------------------------------------------------------------------------------------------------------
  # K1: observed factors, no selection on gm
  # -----------------------------------------------------------------------------------------------------------------
  if (!flag.dgp.issue){
  cat('---------------------------------------- \n')
  cat('K1: known factors, known model \n')
  cat('---------------------------------------- \n \n')
  out_k1=try(fadtwo(y=y,x=x,f=g0, method='joint', L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm=L.gm, U.gm=U.gm, tau1=tau1, tau2=tau2),
             silent = FALSE,
             outFile = getOption("try.outFile", default = stderr())
            )          
  
  if (class(out_k1) != 'try-error') {
  print(out_k1)
  state.k1 = (g0 %*%  out_k1$gm.hat > 0 )
  result.k1[i.rep,] = c( out_k1$bt.hat,  out_k1$dt.hat,  out_k1$gm.hat,  out_k1$objval, state.k1)
  
    # Coverage 1: Use the OLS CI formula and calculate manually
  ap.hat.k1 = c(out_k1$bt.hat, out_k1$dt.hat)
  reg.k1 = cbind(x, x*matrix(state.k1,n.obs,d.x))
  sig2.k1 = sum( (y-reg.k1%*%ap.hat.k1)^2 ) / (n.obs - (2*d.x))
  vcov.k1 = sig2.k1 * solve(t(reg.k1) %*% reg.k1)
  se.ap.k1 = sqrt( diag(vcov.k1) )
  ci1.k1 = cbind(ap.hat.k1 - 1.96*se.ap.k1, ap.hat.k1 + 1.96*se.ap.k1)
  colnames(ci1.k1) = colnames(ci.o)
  cover1.k1[i.rep,] = (ci1.k1[,1] <= c(bt0,dt0)) & (ci1.k1[,2] >= c(bt0,dt0))
  
    # Coverage 2: run the ols again
  ols.k1 = lm(y~reg.k1-1)
  ci2.k1 = confint(ols.k1)
  cover2.k1[i.rep,] = (ci2.k1[,1] <= c(bt0,dt0)) & (ci2.k1[,2] >= c(bt0,dt0))
  } else {
    flag.dgp.issue = TRUE
  }
  }
  # -----------------------------------------------------------------------------------------------------------------
  # K2: known factors, selection on gm
  # -----------------------------------------------------------------------------------------------------------------
  if (!flag.dgp.issue){
  cat('---------------------------------------- \n')
  cat('K2: known factors, unknown model \n')
  cat('---------------------------------------- \n \n')
  g1=g[,1]
  g2=cbind(g[,(2:(feasible.f-1))])
  
  out_k2=try(fadtwo_selection(y=y, x=x, f1=g1, f2=g2, method='joint', L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm1=L.gm1, 
                          U.gm1=U.gm1, L.gm2=L.gm2,U.gm2=U.gm2,L.gm3=L.gm3, U.gm3=U.gm3, L.p=L.p, U.p=U.p, tau1=tau1,tau2=tau2, 
                          eta=eta, params=params, ld=ld),
             silent = FALSE,
             outFile = getOption("try.outFile", default = stderr())
            )
  if (class(out_k2) != 'try-error') {
  print(out_k2)
  state.k2 = (g %*%  out_k2$gm.hat > 0 )
  result.k2[i.rep,] = c(out_k2$bt.hat, out_k2$dt.hat, out_k2$gm.hat, out_k2$objval, state.k2)
  
    # Coverage 1: Use the OLS CI formula and calculate manually
  ap.hat.k2 = c(out_k2$bt.hat, out_k2$dt.hat)
  reg.k2 = cbind(x, x*matrix(state.k2,n.obs,d.x))
  sig2.k2 = sum( (y-reg.k2%*%ap.hat.k2)^2 ) / (n.obs - (2*d.x))
  vcov.k2 = sig2.k2 * solve(t(reg.k2) %*% reg.k2)
  se.ap.k2 = sqrt( diag(vcov.k2) )
  ci1.k2 = cbind(ap.hat.k2 - 1.96*se.ap.k2, ap.hat.k2 + 1.96*se.ap.k2)
  colnames(ci1.k2) = colnames(ci.o)
  cover1.k2[i.rep,] = (ci1.k2[,1] <= c(bt0,dt0)) & (ci1.k2[,2] >= c(bt0,dt0))
  
    # Coverage 2: run the ols again
  ols.k2 = lm(y~reg.k2-1)
  ci2.k2 = confint(ols.k2)
  cover2.k2[i.rep,] = (ci2.k2[,1] <= c(bt0,dt0)) & (ci2.k2[,2] >= c(bt0,dt0))
  } else {
    flag.dgp.issue = TRUE
  }
  }
  
  # -----------------------------------------------------------------------------------------------------------------
  # U: unobserved factors, no model selection, using all factors
  # -----------------------------------------------------------------------------------------------------------------
  if (!flag.dgp.issue){
  cat('---------------------------------------- \n')
  cat('U: unknown factors, all factors included \n')
  cat('---------------------------------------- \n \n')
  pca = get_factors(Xti/sqrt(n.obs*n.crs), feasible.f-1)
  f1 = pca$F
  V_T = diag( pca$eigen.values[c(1:K)] )
  f.hat = cbind(f1[,(1:(K))], -1)
  # Calculate the random gm_0
  g1 = g[,c(1:3)]
  LD = my.dat$Lambda
  H2=t(f1)%*%g1/n.obs
  H3= t(LD)%*%LD/n.crs
  H.bar_T = t(solve(V_T) %*% H2 %*% H3)
  H_T = diag(ncol(H.bar_T)+1)
  H_T[(1:K),(1:K)] = H.bar_T
  gm0 = solve(H_T)%*%phi0

  L.gm.full.u = c(gm0[1], rep(-bnd, K))
  U.gm.full.u = c(gm0[1], rep(bnd, K))
  
  
  out_u = try(fadtwo(y=y,x=x,f=f.hat, method='joint', L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, 
                     L.gm=L.gm.full.u, U.gm=U.gm.full.u, tau1=tau1, tau2=tau2), 
              silent = FALSE,
              outFile = getOption("try.outFile", default = stderr())
              )
  
  if (class(out_u) != 'try-error') {
    print(out_u)
    state.u = (f.hat %*%  out_u$gm.hat > 0 )
    if ( out_u$dt.hat[1] < 0){
      bt.u =  out_u$bt.hat +  out_u$dt.hat
      dt.u = - out_u$dt.hat
      state.u = !(state.u)
    } else {
      bt.u =  out_u$bt.hat
      dt.u =  out_u$dt.hat
    }
    
    result.u[i.rep,] = c(bt.u, dt.u,  out_u$gm.hat,  gm0, out_u$objval, state.u)
    
    
    # Coverage 1: Use the OLS CI formula and calculate manually
    ap.hat.u = c(bt.u, dt.u)
    reg.u = cbind(x, x*matrix(state.u,n.obs,d.x))
    sig2.u = sum( (y-reg.u%*%ap.hat.u)^2 ) / (n.obs - (2*d.x))
    vcov.u = sig2.u * solve(t(reg.u) %*% reg.u)
    se.ap.u = sqrt( diag(vcov.u) )
    ci1.u = cbind(ap.hat.u - 1.96*se.ap.u, ap.hat.u + 1.96*se.ap.u)
    colnames(ci1.u) = colnames(ci.o)
    cover1.u[i.rep,] = (ci1.u[,1] <= c(bt0,dt0)) & (ci1.u[,2] >= c(bt0,dt0))
    
    # Coverage 2: run the ols again
    ols.u = lm(y~reg.u-1)
    ci2.u = confint(ols.u)
    cover2.u[i.rep,] = (ci2.u[,1] <= c(bt0,dt0)) & (ci2.u[,2] >= c(bt0,dt0))  
  } else {
    flag.dgp.issue = TRUE
  }
  }
  
  # If there is an issue in DGP (e.g. positive semidefinite), then we reset the loop number and thraw away the draw.
  if (!flag.dgp.issue) {
    i.rep = i.rep + 1
  }
#  save.image(paste('o_sim',n.crs,'.RData',sep=''))
}





# Model/State Selection Consistency
  # K1
ave.correct.state.k1 = rep(NA, n.rep)
for (i in (1:n.rep)){
  ave.correct.state.k1[i] = mean(result.o[i,paste('state',1:n.obs,sep='')] ==result.k1[i,paste('state',1:n.obs,sep='')])
}
mean.cor.state.k1 = mean(ave.correct.state.k1)
sd.cor.state.k1 = sd(ave.correct.state.k1)

cat('K1: State Prediction =', round(mean.cor.state.k1,3),'(',round(sd.cor.state.k1,3),')','\n')


  # K2
gm.hat.k2 = result.k2[, c('gm1', 'gm2', 'gm3', 'gm4')]

sel.1.k2 = mean( (gm.hat.k2[, 2] == 0) * (gm.hat.k2[, 3] == 0) ) 
sel.2.k2 = mean( (gm.hat.k2[, 2] != 0) * (gm.hat.k2[, 3] == 0) ) 
sel.3.k2 = mean( (gm.hat.k2[, 2] == 0) * (gm.hat.k2[, 3] != 0) ) 
sel.4.k2 = mean( (gm.hat.k2[, 2] != 0) * (gm.hat.k2[, 3] != 0) ) 

ave.correct.state.k2 = rep(NA, n.rep)
for (i in (1:n.rep)){
  ave.correct.state.k2[i] = mean(result.o[i,paste('state',1:n.obs,sep='')] ==result.k2[i,paste('state',1:n.obs,sep='')])
}
mean.cor.state.k2 = mean(ave.correct.state.k2)
sd.cor.state.k2 = sd(ave.correct.state.k2)

cat('K2: Selection (1,4) =', round(sel.1.k2,3),'\n')
cat('K2: Selection (1,2,4) =', round(sel.2.k2,3),'\n')
cat('K2: Selection (1,3,4) =', round(sel.3.k2,3),'\n')
cat('K2: Selection (1,2,3,4) =', round(sel.4.k2,3),'\n')

cat('K2: State Prediction =', round(mean.cor.state.k2,3),'(',round(sd.cor.state.k2,3),')','\n')

  # U
ave.correct.state.u = rep(NA, n.rep)
for (i in (1:n.rep)){
  ave.correct.state.u[i] = mean(result.o[i,paste('state',1:n.obs,sep='')] ==result.u[i,paste('state',1:n.obs,sep='')])
}
mean.cor.state.u = mean(ave.correct.state.u)
sd.cor.state.u = sd(ave.correct.state.u)

cat('K1: State Prediction =', round(mean.cor.state.u,3),'(',round(sd.cor.state.u,3),')','\n')

# Calculate the covarage rate
coverage1.o = apply(cover.o, 2, mean)
coverage1.k1 = apply(cover1.k1, 2, mean)
coverage2.k1 = apply(cover2.k1, 2, mean)
coverage1.k2 = apply(cover1.k2, 2, mean)
coverage2.k2 = apply(cover2.k2, 2, mean)
coverage1.u = apply(cover1.u, 2, mean)
coverage2.u = apply(cover2.u, 2, mean)

cat('------------------ \n')
cat('Coverage Rages \n ')
cat('------------------ \n')
cat('O: Coverage 1  =', coverage1.o, '\n')
cat('K1: Coverage 1 =', coverage1.k1, '\n')
cat('K1: Coverage 2 =', coverage2.k1, '\n')
cat('K2: Coverage 1 =', coverage1.k2, '\n')
cat('K2: Coverage 2 =', coverage2.k2, '\n')
cat('U: Coverage 1  =', coverage1.u, '\n')
cat('U: Coverage 2  =', coverage2.u, '\n')

time.end = proc.time()[3]
runtime = time.end - time.start
cat('\n \n ---------------------------- \n Runtime     =', runtime, 'sec','\n')  



#-------------------------------------------------------------------
# Print out Table 1.1 and 1.2 and save it as '../results/table1.txt'
#-------------------------------------------------------------------

library('xtable')

bt0.m = matrix(1,n.rep,d.x)
dt0.m = matrix(1,n.rep,d.x)
phi0.m = matrix(phi0,n.rep,feasible.f, byrow=T)

# Oracle 1: Know the true state
bt.o = result.o[,paste('bt',(1:d.x),sep='')]
dt.o = result.o[,paste('dt',(1:d.x),sep='')]
mb.bt.o = apply(bt.o - bt0.m, 2, mean)
mb.dt.o = apply(dt.o - dt0.m, 2, mean)
rmse.bt.o = apply(bt.o - bt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.dt.o = apply(dt.o - dt0.m, 2, norm, type="2") / sqrt(n.rep)
re.o = rbind( cbind(mb.bt.o, rmse.bt.o), cbind(mb.dt.o, rmse.dt.o) )
re.o = cbind(re.o, coverage1.o)

# K1: Known factors / Known non-zero gamma (or phi)
bt.k1 = result.k1[,paste('bt',(1:d.x),sep='')]
dt.k1 = result.k1[,paste('dt',(1:d.x),sep='')]
gm.k1 = result.k1[,paste('gm',(1:K),sep='')]
mb.bt.k1 = apply(bt.k1 - bt0.m, 2, mean)
mb.dt.k1 = apply(dt.k1 - dt0.m, 2, mean)
mb.gm.k1 = apply(gm.k1 - phi0.m[,J0], 2, mean)
rmse.bt.k1 = apply(bt.k1 - bt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.dt.k1 = apply(dt.k1 - dt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.gm.k1 = apply(gm.k1 - phi0.m[,J0], 2, norm, type="2") / sqrt(n.rep)
re.k1 = rbind( cbind(mb.bt.k1, rmse.bt.k1), cbind(mb.dt.k1, rmse.dt.k1))
re.k1 = cbind(re.k1, coverage1.k1)
re.k1 = rbind(re.k1, cbind(mb.gm.k1, rmse.gm.k1, NA) )

# K2: Known factors / Unknown non-zero gamma (or phi), i.e. model selection
bt.k2 = result.k2[,paste('bt',(1:d.x),sep='')]
dt.k2 = result.k2[,paste('dt',(1:d.x),sep='')]
gm.k2 = result.k2[,paste('gm',(1:feasible.f),sep='')]
mb.bt.k2 = apply(bt.k2 - bt0.m, 2, mean)
mb.dt.k2 = apply(dt.k2 - dt0.m, 2, mean)
mb.gm.k2 = apply(gm.k2 - phi0.m, 2, mean)
rmse.bt.k2 = apply(bt.k2 - bt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.dt.k2 = apply(dt.k2 - dt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.gm.k2 = apply(gm.k2 - phi0.m, 2, norm, type="2") / sqrt(n.rep)
re.k2 = rbind( cbind(mb.bt.k2, rmse.bt.k2), cbind(mb.dt.k2, rmse.dt.k2) )
re.k2 = cbind(re.k2, coverage1.k2)
re.k2 = rbind(re.k2, cbind(mb.gm.k2, rmse.gm.k2, NA) )

# U: Estimated factors 
bt.u = result.u[,paste('bt',(1:d.x),sep='')]
dt.u = result.u[,paste('dt',(1:d.x),sep='')]
gm.u = result.u[,paste('gm',(1:feasible.f),sep='')]
gm0.u = result.u[,paste('gm0_',(1:feasible.f),sep='')]
gm.u.centered = gm.u - gm0.u
mb.bt.u = apply(bt.u - bt0.m, 2, mean)
mb.dt.u = apply(dt.u - dt0.m, 2, mean)
mb.gm.u = apply(gm.u.centered, 2, mean)
rmse.bt.u = apply(bt.u - bt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.dt.u = apply(dt.u - dt0.m, 2, norm, type="2") / sqrt(n.rep)
rmse.gm.u = apply(gm.u.centered, 2, norm, type="2") / sqrt(n.rep)
re.u = rbind( cbind(mb.bt.u, rmse.bt.u), cbind(mb.dt.u, rmse.dt.u))
re.u = cbind(re.u, coverage1.u)
re.u = rbind(re.u, cbind(mb.gm.u, rmse.gm.u, NA) )

re.tb = round(rbind(re.o, re.k1, re.k2, re.u),4)
colnames(re.tb) = c('Mean Bias', 'RMSE', 'Coverage')


table.1.1 = data.frame(re.tb[c(-9,-16,-24),])
bt.all = c('beta_1', 'beta_2')
dt.all = c('delta_1', 'delta_2')
theta.o = c('theta_2', 'theta_4')
theta.u = c('theta_2', 'theta_3', 'theta_4')
gamma = c('gamma_2', 'gamma_3', 'gamma_4')
par.names = c(bt.all,dt.all, bt.all, dt.all, theta.o, bt.all, dt.all, theta.u, bt.all, dt.all, gamma)
scenarios = c(rep(1,4), rep(2,6), rep(3,7), rep(4,7))
table.1.1 = cbind(scenarios, par.names,table.1.1)

table.1.2 = data.frame(rbind(
  cbind(round(mean.cor.state.k1,4), round(sd.cor.state.k1,4), 'NA'),
  cbind(round(mean.cor.state.k2,4), round(sd.cor.state.k2,4), round(sel.2.k2,4)),
  cbind(round(mean.cor.state.u,4), round(sd.cor.state.u,4), 'NA')
))

colnames(table.1.2) = c('Ave. Cor. Regime Predictions', 's.d.', 'Correct Factor Selection')
scenario2 = c(2,3,4)
table.1.2 = cbind(scenario2, table.1.2)


sink('../results/table1.txt')

print("Table 1 - Part 1 - MB, RMSE, Coverage")
print(table.1.1)
cat('\n', '\n')
print("Table 1 - Part 2 - Regime Prediction")
print(table.1.2)

sink()



save.image(paste('o_sim-table1.RData',sep=''))


