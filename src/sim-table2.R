# File name: sim-table2.R
# 
# A code for simulation studies in "Factor-driven Two-regime Regression"
# : This code generates Table 2 in the paper and save it as "../results/table2.txt"
#
# Last Update
#   2018-09-26
#

rm(list=ls())
time.start = proc.time()[3]
n.sim='001'
set.seed(n.sim)

source('lib_fadtwo.R')

dgp3 <- function(n.obs,n.crs, d.x, rho.g, rho.x, K, sg_eps, bt, dt, phi){

  # Generating the Factor processes
    # g_j,t = rho.g_j * g_j,t-1 + u_j,t  j = 1,...,K
    u = matrix(rnorm(K*(n.obs+1)), K, n.obs+1)
    g1 = matrix(NA, K, n.obs)
    g1 = cbind(rep(0,K), g1)    # Add the inital value
    n.rho.g = length(rho.g)
    m.rho.g = matrix(0, n.rho.g, n.rho.g)
    diag(m.rho.g) = rho.g
    for (i in c(2:(n.obs+1))){
        g1[,i] = m.rho.g %*% g1[,i-1] + u[,i]  
    }
    g1 = g1[,c(2:(n.obs+1)), drop=FALSE]

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

  return(list(y=y, x=x, g=g, state0=state))  
}

generate_Yti = function(g1, n.crs, rho.e){
  K = ncol(g1)
  n.obs = nrow(g1)
  
  # Generating the Factor processes
  # (1) e_it = rho_i * e_it-1 + eps_it
  eps_e = matrix(rnorm(n.crs*(n.obs+1)), n.crs, n.obs+1)
  e = matrix(NA, n.crs, n.obs)
  e = cbind(rep(0,n.crs), e)    # Add the inital value
  for (i in c(2:(n.obs+1))){
    e[,i] = diag(rho.e) %*% e[,i-1] + eps_e[,i]
  }
  e = e[,c(2:(n.obs+1))]
  # (3) X = Lambda * g1 + \sqrt{K}*e
  Lambda = matrix(rnorm(n.crs*K)*sqrt(K), n.crs, K)
  Yit = Lambda %*% t(g1) + sqrt(K) * e
  
  return(list(Yti=t(Yit), Lambda=Lambda))
}

# parameters for DGP
n.obs = 200     
crs.set = c(100, 200, 400, 1600)
n.crs.set = length(crs.set)
n.rep = 1000
result.u = length(crs.set)


d.x = 2             # Dimension of x_t
rho_x = 0.5         # AR(1) coefficient for x_t
K = 1               # The number of non-constant factors
sg_eps = 0.5        # Standard deviation of epsilon (the error term of the model)

p = 0
feasible.f = 2+p  # This dimension includes the constant term (f1,f2,f3,-1)

rho_e = list()
for (i.crs in (1:n.crs.set)){
  n.crs = crs.set[i.crs]
  rho_e[[i.crs]] = runif(n.crs, min=0.3, max=0.5)
}
rho_g = runif(K, min=0.2, max=0.8)

bt0 = rep(1,d.x)
dt0 = rep(1,d.x)
ap0 = c(bt0, dt0)
phi0 = c(1,2/3)


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






# Deinfe a list object to save the results: for all crs.set
result.u = list()
cover.u = list()
state0.u = list()

coverage.u = list()
ave.correct.state.u = list()
mean.cor.state.u = list()
sd.cor.state.u = list()

for (i.crs in (1:n.crs.set)){
  result.u[[i.crs]] = matrix(NA, n.rep, d.x+d.x+2*feasible.f+1+1+n.obs)
  colnames(result.u[[i.crs]]) = c(paste('bt',(1:d.x),sep=''), paste('dt',(1:d.x),sep=''), paste('gm',(1:feasible.f),sep=''), paste('gm0_',(1:feasible.f),sep=''), 'gm.ratio','objval', paste('state',(1:n.obs),sep=''))
  cover.u[[i.crs]] = matrix(NA, n.rep, d.x*2)
  ave.correct.state.u[[i.crs]] = rep(NA, n.rep)
}


i.rep=1


while (i.rep <= n.rep) {
  flag.dgp.issue = FALSE
  cat('----------------------- \n', i.rep, 'iterations','\n','-------------------------','\n')
  # Data Generation for y, x, and g
  my.dat = dgp3(n.obs,n.crs, d.x, rho_g, rho_x, K, sg_eps, bt=bt0, dt=dt0, phi=phi0) 
  
  y=my.dat$y
  x=my.dat$x
  g=my.dat$g
  g1=g[,(1:K), drop=FALSE]
  state0=my.dat$state0
  J0 = (phi0 != 0)
  g0=(my.dat$g)[, J0]
  
  
  # -----------------------------------------------------------------------------------------------------------------
  # U: unobserved factors, no model selection, using all factors
  # -----------------------------------------------------------------------------------------------------------------
  cat('---------------------------------------- \n')
  cat('U: unknown factors, all factors included \n')
  cat('---------------------------------------- \n \n')
  
  for (i.crs in (1:n.crs.set)){
    # Generate Y_ti for unobserved factors
    n.crs = crs.set[i.crs]
    gen_Yti = generate_Yti(g1, n.crs, rho_e[[i.crs]])
    Y.ti = gen_Yti$Yti
    LD = gen_Yti$Lambda
    cat('------------------------------------------------------------- \n')
    cat('N =',n.crs,'\n')
    cat('------------------------------------------------------------- \n')

    if (!flag.dgp.issue){
      pca = get_factors(Y.ti/sqrt(n.obs*n.crs), feasible.f-1)
      f1 = as.matrix(pca$F)
      n.eigen=length(pca$eigen.values[c(1:K)])
      V_T = diag(n.eigen)
      diag(V_T) = pca$eigen.values[c(1:K)]
      f.hat = cbind(f1[,(1:(K))], -1)
      # Calculate the random gm_0
      H2=t(f1)%*%g1/n.obs
      H3= t(LD)%*%LD/n.crs
      H.bar_T = t(solve(V_T) %*% H2 %*% H3)
      H_T = diag(ncol(H.bar_T)+1)
      H_T[(1:K),(1:K)] = H.bar_T
      gm0 = solve(H_T)%*%phi0
      
      L.gm.full.u = c(L.gm1, rep(-bnd, K))
      U.gm.full.u = c(U.gm1, rep(bnd, K))
      
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
        
        result.u[[i.crs]][i.rep,] = c(bt.u, dt.u, out_u$gm.hat,  gm0, gm0[2,1]/gm0[1,1], out_u$objval, state.u)
        
        # Calculate the average correct state prediction rates
        ave.correct.state.u[[i.crs]][i.rep] = mean(state0 ==result.u[[i.crs]][i.rep,paste('state',1:n.obs,sep='')])
        
        # Coverage 1: Use the OLS CI formula and calculate manually
        ap.hat.u = c(bt.u, dt.u)
        reg.u = cbind(x, x*matrix(state.u,n.obs,d.x))
        sig2.u = sum( (y-reg.u%*%ap.hat.u)^2 ) / (n.obs - (2*d.x))
        vcov.u = sig2.u * solve(t(reg.u) %*% reg.u)
        se.ap.u = sqrt( diag(vcov.u) )
        ci1.u = cbind(ap.hat.u - 1.96*se.ap.u, ap.hat.u + 1.96*se.ap.u)
        cover.u[[i.crs]][i.rep,] = (ci1.u[,1] <= c(bt0,dt0)) & (ci1.u[,2] >= c(bt0,dt0))
      } else {
        flag.dgp.issue = TRUE
        error.dgp = cbind(y, x, g)
        error.code = out_u
        error.yti = Y.ti
      }
    }
    
  }

  
    
    
  # If there is an issue in DGP (e.g. positive semidefinite), then we reset the loop number and thraw away the draw.
  if (!flag.dgp.issue) {
    i.rep = i.rep + 1
  }
  save.image(paste('o_sim',n.crs,'.RData',sep=''))
}


# Calculate the covarage rate

for (i.crs in (1:n.crs.set)){
  coverage.u[[i.crs]] = apply(cover.u[[i.crs]], 2, mean)
  
  #for (i in (1:n.rep)){
  #  ave.correct.state.u[[i.crs]][i] = mean(state0 ==result.u[[i.crs]][i,paste('state',1:n.obs,sep='')])
  #}
  mean.cor.state.u[[i.crs]] = mean(ave.correct.state.u[[i.crs]] )
  sd.cor.state.u[[i.crs]] = sd(ave.correct.state.u[[i.crs]] )
  cat('N=',crs.set[i.crs],'State Prediction =', round(mean.cor.state.u[[i.crs]],3),'(',round(sd.cor.state.u[[i.crs]],3),')','\n')
  
}


time.end = proc.time()[3]
runtime = time.end - time.start
cat('\n \n ---------------------------- \n Runtime     =', runtime, 'sec','\n')  

#----------------------------------------------------------
# Print out Table 2 and save it as '../results/table2.txt' 
#----------------------------------------------------------

bt0.m = matrix(1,n.rep,d.x)
dt0.m = matrix(1,n.rep,d.x)

sink('../results/table2.txt')
cat('---------------------------------------','\n')

# U: Estimated factors 
for (i.crs in (1:n.crs.set)){
  bt.u = result.u[[i.crs]][,paste('bt',(1:d.x),sep='')]
  dt.u = result.u[[i.crs]][,paste('dt',(1:d.x),sep='')]
  gm.u = result.u[[i.crs]][,paste('gm',(1:feasible.f),sep='')]
  gm.ratio = result.u[[i.crs]][,'gm.ratio']
  gm.u.ratio = gm.u[,2] - gm.ratio
  mb.bt.u = apply(bt.u - bt0.m, 2, mean)
  mb.dt.u = apply(dt.u - dt0.m, 2, mean)
  mb.gm.u = mean(gm.u.ratio)
  rmse.bt.u = apply(bt.u - bt0.m, 2, norm, type="2") / sqrt(n.rep)
  rmse.dt.u = apply(dt.u - dt0.m, 2, norm, type="2") / sqrt(n.rep)
  rmse.gm.u = norm(gm.u.ratio, type="2") / sqrt(n.rep)
  re.u = rbind( cbind(mb.bt.u, rmse.bt.u), cbind(mb.dt.u, rmse.dt.u))
  re.u = cbind(re.u, coverage.u[[i.crs]])
  re.u = rbind(re.u, cbind(mb.gm.u, rmse.gm.u, NA) )
  
  re.tb = rbind(re.u)

  table.2 = data.frame(round(re.tb[,c(1:2)],4))
  table.2 = cbind(c('beta_1','beta_2','delta_1','delta_2','gamma2/gamma1'), table.2)
  colnames(table.2) = c('Par', 'Mean Bias', 'RMSE')
  rownames(table.2) = NULL

  cat('N = ', crs.set[i.crs],'\n')
  print(table.2)
  cat('---------------------------------------','\n')
  cat(cbind('Ave. Prediction   ',round(mean.cor.state.u[[i.crs]],4), '(', round(sd.cor.state.u[[i.crs]],4)),')','\n')
  cat('---------------------------------------','\n')
}

sink()




save.image(paste('../results/o-sim-table2.RData',sep=''))






