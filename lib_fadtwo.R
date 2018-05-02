# File Name: lib_FADTWO.R
# This library contains all functions necessary for "FActor Driven TWO-regime regresssion (FADTWO)"
#
# Latest Update: 
#   2018-02-22	  Merge with model selection codes
#   2018-01-27    Original Code
# 


#------------------------------------------------------------------------------------------------------------------------------------
# Function name: fadtwo
#
# This function estimates the model WITHOUT the model selection penalty term
# 

fadtwo <- function(y,x,f, method="joint", L.bt=NULL, U.bt=NULL, L.dt=NULL, U.dt=NULL, L.gm=NULL, U.gm=NULL, tau1, tau2, 
                   eta=1e-6, params=list(OutputFlag=1, FeasibilityTol=1e-9),
                   grid=NULL, max.iter=2) {
  
  # Model: y = x'bt + x'dt * 1{f'gm > 0} + eps
  #
  # Input:
  #   y: outcome variable
  #   x: covariates
  #   f: factors
  #   L.bt / U.bt: Lower and upper bounds for bt. The dim should be equal to the dim of bt
  #   L.dt / U.dt: Lower and upper bounds for dt
  #   L.gm / U.gm: Lower and upper bounds for gm
  #   tau1: Lower bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   tau2: Upper bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   eta: effective zero
  #   params: parameters for gurobi engine
  #
  x = as.matrix(x)
  f = as.matrix(f)
  n.obs = nrow(x)     # Number of observations
  d.x = ncol(x)       # Dimension of regressors
  d.f = ncol(f)       # Dimension of factors
  
  if (method=="joint"){
    # Calculate Big-M  
    A.gm = c(1,-1) %x% diag(rep(1,d.f)) 
    b.gm = c(U.gm, -L.gm)
    M = rep(NA,n.obs)
    for (i in (1:n.obs)){
      M[i] = get_m(f=f[i,],A=A.gm,b=b.gm)$m
    }  
    
    # Estimate parameters by solving the mixed integer quadratic programming
    # 1. Construct the objective function  
    L.obj = get_L(x=x,y=y,L.dt=L.dt, d.f=d.f)
    Q.obj = get_Q(x=x, L.dt=L.dt, d.f=d.f)
    objcon = mean(y^2)
    
    # 2. Build the constraints
    const = build_constraint(L.bt, U.bt, L.dt, U.dt, L.gm, U.gm, M=M, eta=eta, d.x=d.x, d.f=d.f, n.obs=n.obs, f=f, tau1=tau1, tau2=tau2)
    A.const = const$A.const
    b.const = const$b.const
    
    # 3. Estimate the model
    result = estimate(y=y, x=x, f=f, Q.obj=Q.obj, L.obj=L.obj, objcon=objcon, A.const=A.const, b.const=b.const, L.bt=L.bt, L.gm=L.gm, params=params)
    opt.par = result$x
    bt.hat = opt.par[1:d.x]
    l.hat = opt.par[(d.x+1):(d.x*(n.obs+1))]
    d.hat = opt.par[(d.x*(n.obs+1) +1):(d.x*(n.obs+1) + n.obs)]
    dt.tilde = opt.par[((d.x*(n.obs+1))+ n.obs + 1):( (d.x*(n.obs+1))+ n.obs  + d.x )]
    # Note that dt.tilde = dt.hat - L.dt.  
    dt.hat = dt.tilde + L.dt  
    gm.hat = opt.par[((d.x*(n.obs+1))+n.obs  + d.x + 1):((d.x*(n.obs+1))+n.obs  + d.x  + d.f)]
    objval=get_objval(y,x,f,bt.hat,dt.hat,gm.hat, eta=eta)
    cat('-----------------------------------------', '\n')
    cat('Results', '\n')
    cat('-----------------------------------------', '\n')
    cat('Estimation Method =',method, '\n')
    cat('Sample size =', n.obs, '\n')
    cat('Obj val     =', objval, '\n')
    cat('Beta.hat   =', bt.hat, '\n')
    cat('Delta.hat   =', dt.hat, '\n')
    cat('Gamma.hat   =', gm.hat, '\n')
    cat('-----------------------------------------', '\n')
    # return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, result.out=result))  # Output for detailed results
    return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, objval=objval))
  }
  else if(method == "iter"){
    # Grid Search for gm
    n.grid=nrow(grid)
    step1.out = step1_grid(y=y,x=x,f=f,grid=grid, eta=eta)
    ap.hat.step1 = step1.out$ap.hat
    gm.hat.step1 = step1.out$gm.hat
    cat('-----------------------------------------', '\n')
    cat('Number of Grid Points = ', n.grid,'\n')
    cat('ap.hat = ', ap.hat.step1,'\n')
    cat('gm.hat = ', gm.hat.step1,'\n')
    cat('-----------------------------------------', '\n')
    
    # Calculate Big-M  
    A.gm = c(1,-1) %x% diag(rep(1,d.f)) 
    b.gm = c(U.gm, -L.gm)
    M = rep(NA,n.obs)
    for (i in (1:n.obs)){
      M[i] = get_m(f=f[i,],A=A.gm,b=b.gm)$m
    } 
    
    # Estimate gm.hat by MIO and Update ap.hat
    bt.pre = ap.hat.step1[c(1:d.x)]
    dt.pre = ap.hat.step1[-c(1:d.x)]
    for (cnt.it in (1:max.iter)){
      step2_1.out = estimate_gm(y=y, x=x, f=f, bt=bt.pre, dt=dt.pre, A=A.gm, b=b.gm, 
                              M=M, tau1=tau1, tau2=tau2, params=params, eta = eta)
    gm.hat = step2_1.out$gm
    d.hat = step2_1.out$d.t
    
    # Update ap.hat
    step2_2.out = estimate_bt_dt(y=y, x=x, f=f, gm=gm.hat, eta=eta)
    bt.hat = step2_2.out$bt.hat
    dt.hat = step2_2.out$dt.hat
    objval=get_objval(y,x,f,bt.hat,dt.hat,gm.hat, eta=eta)
    cat('-----------------------------------------', '\n')
    cat('Results', '\n')
    cat('-----------------------------------------', '\n')
    cat('Estimation Method =',method, '\n')
    cat('Number of Iterations = ', cnt.it, '\n')
    cat('Sample size =', n.obs, '\n')
    cat('Obj val     =', objval, '\n')
    cat('Beta.hat    =', bt.hat, '\n')
    cat('Delta.hat   =', dt.hat, '\n')
    cat('Gamma.hat   =', gm.hat, '\n')
    cat('-----------------------------------------', '\n')
    #return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, result.out=step2_1.out)) # Output for detailed results
    }
    bt.pre=bt.hat
    dt.pre=dt.hat
    return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, objval=objval))
  }
  else {
    stop("The 'method' option should be either 'joint' or 'iter'.")
  }
  
  
  
  
}

#------------------------------------------------------------------------------------------------------------
# Function: fadtwo_selection
# 
# Estimate the model when the objective function includes the \ell_0 penalty for model selection
#
#

fadtwo_selection <- function(y, x, f1, f2, method, L.bt=NULL, U.bt=NULL, L.dt=NULL, U.dt=NULL, L.gm1, U.gm1, L.gm2, U.gm2, L.gm3, U.gm3, 
                             L.p, U.p, tau1, tau2, eta, params=list(OutputFlag=1, FeasibilityTol=1e-9), grid=NULL, max.iter=2, ld) 


  {

  # Model: y = x'bt + x'dt * 1{f'gm > 0} + eps
  #
  # Input:
  #   y: outcome variable
  #   x: covariates
  #   f1: factors known to be active 
  #   f2: factors to be selected
  #   method: estimation method, "joint" or "iter"
  #   L.bt / U.bt: Lower and upper bounds for bt. The dim should be equal to the dim of bt
  #   L.dt / U.dt: Lower and upper bounds for dt
  #   L.gm1 / U.gm1: Lower and upper bounds for gm1 (corresponding to f1)
  #   L.gm2 / U.gm2: Lower and upper bounds for gm2 (corresponding to f2)
  #   L.gm3 / U.gm3: scalars, Lower and upper bound for gm3 (parameter for the constant term -1)
  #   L.p: the minimum number of factors to be included among f2 ('p lower bar' in the paper)
  #   U.p: the maximum number of factors to be included among f2 ('p upper bar' in the paper)
  #   tau1: Lower bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   tau2: Upper bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   eta: effective zero
  #   params: parameters for gurobi engine
  #   grid: grid points for searching gm.hat when method="iter"
  #   max.iter: the number of iterations (updates) after the grid search when method="iter"
  #   ld: a constant multiplied to the penalty term, AIC/BIC coefficient ('lambda' in the paper)
  #

  x = as.matrix(x)
  f1 = as.matrix(f1)
  f2 = as.matrix(f2)
  f = cbind(f1, f2, -1)
  n.obs = nrow(x)     # Number of observations
  d.x = ncol(x)       # Dimension of x_t
  d.f = ncol(f)       # Dimension of f_t
  p = ncol(f2)			 # Dimension of f_2t
  
  if (method=="joint"){
    # Calculate Big-M  
    A.gm = c(1,-1) %x% diag(rep(1,d.f)) 
    b.gm = c(U.gm1, U.gm2, U.gm3, -L.gm1, -L.gm2, -L.gm3)
    M = rep(NA,n.obs)
    for (i in (1:n.obs)){
      M[i] = get_m(f=f[i,],A=A.gm,b=b.gm)$m
    }  
    
    # Estimate parameters by solving the mixed integer quadratic programming
    # 1. Construct the objective function  
    L.obj = get_L_selection(x=x,y=y,L.dt=L.dt, d.f=d.f, ld=ld, p=p)
    Q.obj = get_Q_selection(x=x, L.dt=L.dt, d.f=d.f, p=p)
    objcon = mean(y^2)
    
    # 2. Build the constraints
    const = build_constraint_selection(L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm1=L.gm1, U.gm1=U.gm1, L.gm2=L.gm2, U.gm2=U.gm2, L.gm3=L.gm3, U.gm3=U.gm3, L.p=L.p, U.p=U.p, M=M, eta=eta, 
                             d.x=d.x, d.f=d.f, n.obs=n.obs, f=f, tau1=tau1, tau2=tau2, p=p)
    A.const = const$A.const
    b.const = const$b.const
    
    # 3. Estimate the model
    result = estimate_selection(y=y, x=x, f=f, Q.obj=Q.obj, L.obj=L.obj, objcon=objcon, A.const=A.const, b.const=b.const, L.bt=L.bt, L.gm=L.gm, params=params, p=p)
    opt.par = result$x
	   names(opt.par) = c(paste('bt',c(1:d.x),sep='_'), paste('l',c(1:(d.x*n.obs)),sep='_'), paste('d',c(1:n.obs),sep='_'), 
								paste('dt.tilde',c(1:d.x),sep='_'), paste('gm',c(1:d.f),sep='_'), paste('e',c(1:p),sep='_'))
    bt.hat = opt.par[paste('bt',c(1:d.x),sep='_')]
    l.hat = opt.par[paste('l',c(1:(d.x*n.obs)),sep='_')]
    d.hat = opt.par[paste('d',c(1:n.obs),sep='_')]
    dt.tilde = opt.par[paste('dt.tilde',c(1:d.x),sep='_')]
    dt.hat = dt.tilde + L.dt            # Note that dt.tilde = dt.hat - L.dt.  
    gm.hat = opt.par[paste('gm',c(1:d.f),sep='_')]
    gm.hat = gm.hat * (gm.hat > eta)    # Only keep gm.hat bigger than the effective zero.
      gm1.hat = gm.hat[c(1:ncol(f1))]
      gm2.hat = gm.hat[c((ncol(f1)+1):(length(gm.hat)-1))]
      gm3.hat = gm.hat[length(gm.hat)]
    e.hat = opt.par[paste('e',c(1:p),sep='_')]
    #objval=get_objval(y,x,f,bt.hat,dt.hat,gm.hat, eta=eta)
    objval=result$objval
    cat('-----------------------------------------', '\n')
    cat('Results', '\n')
    cat('-----------------------------------------', '\n')
    cat('Estimation Method =',method, '\n')
    cat('Sample size =', n.obs, '\n')
    cat('Obj val     =', objval, '\n')
    cat('Beta.hat   =', bt.hat, '\n')
    cat('Delta.hat   =', dt.hat, '\n')
    cat('Gamma.hat   =', gm.hat, '\n')
    cat(' Gamma1.hat =', gm1.hat, '\n')
    cat(' Gamma2.hat =', gm2.hat, '\n')
    cat(' Gamma3.hat =', gm3.hat, '\n')
    cat('-----------------------------------------', '\n')
    # return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, result.out=result))  # Output for detailed results
    return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, objval=objval))
  }
  else if(method == "iter"){
    # Grid Search for gm
    n.grid=nrow(grid)
    step1.out = step1_grid(y=y,x=x,f=f,grid=grid, eta=eta)
    ap.hat.step1 = step1.out$ap.hat
    gm.hat.step1 = step1.out$gm.hat
    cat('-----------------------------------------', '\n')
    cat('Number of Grid Points = ', n.grid,'\n')
    cat('ap.hat = ', ap.hat.step1,'\n')
    cat('gm.hat = ', gm.hat.step1,'\n')
    cat('-----------------------------------------', '\n')

    # Calculate Big-M  
    A.gm = c(1,-1) %x% diag(rep(1,d.f)) 
    b.gm = c(U.gm1, U.gm2, U.gm3, -L.gm1, -L.gm2, -L.gm3)
    M = rep(NA,n.obs)
    for (i in (1:n.obs)){
      M[i] = get_m(f=f[i,],A=A.gm,b=b.gm)$m
    } 
    
    bt.pre = ap.hat.step1[c(1:d.x)]
    dt.pre = ap.hat.step1[-c(1:d.x)]
    
    # Test to compare it with joint estimation codes
    # bt.pre=c(0.1137517, 0.1308253, 0.09609715, -0.02124596, -0.009057051, 0.01754676, 1.407947, -0.9421353, 0.7915665, -0.3855972, 0.4666518, -0.7873762, 0.3730487, 0.02099868)
    # dt.pre=c(-0.09237944, -0.112688, -0.1321168, 0.03301762, -0.03461472, 0.02713847, -0.124312, 0.6114814, -0.5539258, 0.1700254, 0.2375884, 0.3583261, -0.9860771, 0.3360492)


    for (cnt.it in (1:max.iter)){
      # Estimate gm.hat by MIO and Update ap.hat
      step2_1.out = estimate_gm_selection(y=y, x=x, f1=f1, f2=f2, bt=bt.pre, dt=dt.pre, L.gm1=L.gm1, U.gm1=U.gm1, L.gm2=L.gm2, 
                                          U.gm2=U.gm2, L.gm3=L.gm3, U.gm3=U.gm3, M=M, tau1=tau1, tau2=tau2, params=params, eta = eta, ld=ld, p=p)
      gm.hat = step2_1.out$gm
      gm.hat = gm.hat * (gm.hat > eta)
      gm1.hat = gm.hat[c(1:ncol(f1))]
		  gm2.hat = gm.hat[c((ncol(f1)+1):(length(gm.hat)-1))]
		  gm3.hat = gm.hat[length(gm.hat)]
      e.hat = step2_1.out$e
      d.hat = step2_1.out$d.t
      
      # Update ap.hat
      step2_2.out = estimate_bt_dt(y=y, x=x, f=f, gm=gm.hat, eta=eta)
      bt.hat = step2_2.out$bt.hat
      dt.hat = step2_2.out$dt.hat
      objval=get_objval(y,x,f,bt.hat,dt.hat,gm.hat, eta=eta) + ld*sum(e.hat)
      cat('-----------------------------------------', '\n')
      cat('Results', '\n')
      cat('-----------------------------------------', '\n')
      cat('Estimation Method =',method, '\n')
      cat('Number of Iterations = ', cnt.it, '\n')
      cat('Sample size =', n.obs, '\n')
      cat('Obj val     =', objval, '\n')
      cat('Beta.hat    =', bt.hat, '\n')
      cat('Delta.hat   =', dt.hat, '\n')
      cat('Gamma.hat   =', gm.hat, '\n')
      cat(' Gamma1.hat =', gm1.hat, '\n')
      cat(' Gamma2.hat =', gm2.hat, '\n')
      cat(' Gamma3.hat =', gm3.hat, '\n')
      cat('-----------------------------------------', '\n')
      
      # Update the initial values of bt.hat and dt.hat
      bt.pre = bt.hat
      dt.pre = dt.hat
    }
    #return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, result.out=step2_1.out)) # Output for detailed results
    return(list(bt.hat=bt.hat, dt.hat=dt.hat, gm.hat=gm.hat, objval=objval))
  }
  else {
    stop("The 'method' option should be either 'joint' or 'iter'.")
  }
  

  

}



# Function Name: get_m
#
# This function constructs 'big-M' which is used for bounds in MIO
#
get_m <-function(f,A,b,flag.print=FALSE){
  
  # Call necessary libraries
  library('gurobi')   # Gurobi library for Linear Programming
  
  # Declare parameters
  dim.gm = length(as.numeric(f))    # The dimension of gamma
  dim.b = length(b)                 # The number of constraint for parameter space
  params <- list(OutputFlag=0)      # Parameters for Gurobi
  
  # Declare a model for positive side
  model.pos <- list()
  
  A.1 = rbind( A, -f )
  model.pos$A          = A.1
  model.pos$obj        = f
  model.pos$modelsense = "max"
  model.pos$rhs        = c(b,0)
  model.pos$sense      = c(rep('<=',dim.b+1))
  model.pos$lb         = c(rep(-10^5,dim.gm))   # w/o lb Gurobi automatically impose 0 as a lower bound. 
  
  result.pos <- gurobi(model.pos, params)
  
  # Declare a model for negative side
  model.neg <- list()
  
  A.2 = rbind( A, f )
  model.neg$A          = A.2
  model.neg$obj        = -f
  model.neg$modelsense = "max"
  model.neg$rhs        = c( b , 0)
  model.neg$sense      = c(rep('<=',dim.b+1))
  model.neg$lb         = c(rep(-10^5,dim.gm))    # Lower bound for gamma is set to be -10^10 to allow negative values. 
  
  result.neg <- gurobi(model.neg, params)
  
  if (flag.print){
    cat('---------------------','\n')
    cat('Solutions','\n')
    cat('---------------------','\n')
    cat('objective value (positive):', result.pos$objval,'\n')
    cat('maximizer (positive):', result.pos$x,'\n')
    cat('objective value (negative):', result.neg$objval,'\n')
    cat('maximizer (negative):', result.neg$x,'\n')
  }
  
  # Check if the solution exists; then relocate the value by add/subtract the constant value f[1]
  if (!is.numeric(result.pos$objval)){
    result.pos$objval = 0
  }
  if (!is.numeric(result.neg$objval)){
    result.neg$objval = 0
  }
  
  if ( result.neg$objval > result.pos$objval ) {
    result = result.neg
  } else {
    result = result.pos
  }
  
  # Raise a flag when the result is not optmal
  if (result$status!="OPTIMAL") {
    cat('Warning: the result is not optimal!','\n')
  }
  
  return(list(m=result$objval, sol=result$x))
}


# Function Name: get_Q
#
# This function constructs a square matrix 'Q' for the quadratic term of the objective function
#
get_Q <-function(x, L.dt, d.f){
  library('Matrix')
  
  T = nrow(as.matrix(x))
  d.x = ncol(as.matrix(x))
  
  outer.x = array(NA, dim=c(d.x,d.x,T)) # d.x x d.x x T arrays to store column outer products
  
  for (i in c(1:T)){
    outer.x[,,i] = x[i,] %*% t(x[i,])
  }
  
  x11 = outer.x[,,1]
  Q.row = x11
  Q.row.L = x11 %*% L.dt
  Q.diag = x11
  Q.diag.L = x11 %*% L.dt
  Q.diag.LL = t(L.dt) %*% x11 %*% L.dt
  for (i in c(2:T)){
    Q.row = cbind(Q.row, outer.x[,,i]) 
    Q.row.L = cbind(Q.row.L, outer.x[,,i]%*%L.dt)
    Q.diag = bdiag(Q.diag,outer.x[,,i])
    Q.diag.L = bdiag(Q.diag.L,outer.x[,,i] %*% L.dt)
    Q.diag.LL = bdiag(Q.diag.LL,t(L.dt)%*%outer.x[,,i] %*% L.dt)
  }
  
  Q1 = cbind(t(x)%*%x, Q.row, Q.row.L)   # The first d.x x (T+1)d.x matrix
  Q2 = cbind(t(Q.row), Q.diag, Q.diag.L)  # Next rows of Q
  Q3 = cbind(t(Q.row.L), t(Q.diag.L), Q.diag.LL)  # Next rows of Q
  Q = rbind(Q1,Q2, Q3) / T
  
  # Add additional 0 columns
  Q = cbind(Q, matrix(0, nrow=nrow(Q), ncol=(d.x+d.f)))
  # Add additional 0 rows
  Q = rbind(Q, matrix(0, nrow=(d.x+d.f), ncol=ncol(Q)))
  
  return(Q)
}

# Function Name: get_Q_selection
#
# This function constructs a square matrix 'Q' for the quadratic term of the objective function with model selection
#
get_Q_selection <-function(x, L.dt, d.f, p){
  library('Matrix')
  
  T = nrow(as.matrix(x))
  d.x = ncol(as.matrix(x))
  
  outer.x = array(NA, dim=c(d.x,d.x,T)) # d.x x d.x x T arrays to store column outer products
  
  for (i in c(1:T)){
    outer.x[,,i] = x[i,] %*% t(x[i,])
  }
  
  x11 = outer.x[,,1]
  Q.row = x11
  Q.row.L = x11 %*% L.dt
  Q.diag = x11
  Q.diag.L = x11 %*% L.dt
  Q.diag.LL = t(L.dt) %*% x11 %*% L.dt
  for (i in c(2:T)){
    Q.row = cbind(Q.row, outer.x[,,i]) 
    Q.row.L = cbind(Q.row.L, outer.x[,,i]%*%L.dt)
    Q.diag = bdiag(Q.diag,outer.x[,,i])
    Q.diag.L = bdiag(Q.diag.L,outer.x[,,i] %*% L.dt)
    Q.diag.LL = bdiag(Q.diag.LL,t(L.dt)%*%outer.x[,,i] %*% L.dt)
  }
  
  Q1 = cbind(t(x)%*%x, Q.row, Q.row.L)   # The first d.x x (T+1)d.x matrix
  Q2 = cbind(t(Q.row), Q.diag, Q.diag.L)  # Next rows of Q
  Q3 = cbind(t(Q.row.L), t(Q.diag.L), Q.diag.LL)  # Next rows of Q
  Q = rbind(Q1,Q2, Q3) / T
  
  # Add additional 0 columns
  Q = cbind(Q, matrix(0, nrow=nrow(Q), ncol=(d.x+d.f+p)))
  # Add additional 0 rows
  Q = rbind(Q, matrix(0, nrow=(d.x+d.f+p), ncol=ncol(Q)))
  
  return(Q)
}


# Function Name: get_L
#
# This function constructs a vector 'L' for the linear term of the objective function
#
get_L <-function(x, y, L.dt, d.f){
  
  n.obs = length(y)
  y = as.matrix(y,n.obs,1)
  d.x = ncol(as.matrix(x))
  
  A.1 = t(y)%*%x
  A.2 = y[1]%*%x[1,]
  A.3 = y[1]*(x[1,]%*%L.dt)
  for (i in c(2:n.obs)){
    A.2 = cbind(A.2, y[i]%*%x[i,])
    A.3 = cbind(A.3, y[i]*(x[i,]%*%L.dt))
  }
  
  A = cbind(A.1,A.2,A.3)
  A = (-2/n.obs) * A
  L = cbind(A, matrix(0,1,(d.x+d.f)))
  
  return(L)
}

# Function Name: get_L_selection
#
# This function constructs a vector 'L' for the linear term of the objective function with model selection
#
get_L_selection <-function(x, y, L.dt, d.f, ld, p){
  
  n.obs = length(y)
  y = as.matrix(y,n.obs,1)
  d.x = ncol(as.matrix(x))
  
  A.1 = t(y)%*%x
  A.2 = y[1]%*%x[1,]
  A.3 = y[1]*(x[1,]%*%L.dt)
  for (i in c(2:n.obs)){
    A.2 = cbind(A.2, y[i]%*%x[i,])
    A.3 = cbind(A.3, y[i]*(x[i,]%*%L.dt))
  }
  
  A = cbind(A.1,A.2,A.3)
  A = (-2/n.obs) * A
  L = cbind(A, matrix(0,1,(d.x+d.f)), matrix(ld,1,p))
  
  return(L)
}



# Function Name: build_constraint
#
# This function constructs the whole linear constraints: Matrix of the LHS and a vector of RHS will be generated
#
build_constraint = function(L.bt, U.bt, L.dt, U.dt, L.gm, U.gm, M, eta, d.x, d.f, n.obs, f, tau1, tau2){
  
  # Build the RHS vector of constatins
  
  # constraints #1
  b.const = c(U.bt, -L.bt)
  # constraints #2
  b.const = c(b.const, U.gm, -L.gm)
  # constraints #3
  b.const = c(b.const, (U.dt - L.dt), rep(0,d.x))
  # constraints #4
  b.const = c(b.const, rep(0,d.x*n.obs))
  # constraints #4.5
  b.const = c(b.const, rep(0,d.x*n.obs))
  # constraints #5
  b.const = c(b.const, M+eta)
  # constraints #6
  b.const = c(b.const, rep(0,n.obs))
  # constraints #7
  b.const = c(b.const, rep(0,n.obs))
  # constraints #7.5
  b.const = c(b.const, rep(0,n.obs))
  # constraints #8
  b.const = c(b.const, rep(sum(U.dt - L.dt),n.obs))
  # constraints #8.5
  b.const = c(b.const, rep(0,n.obs))
  # constraints #9
  b.const = c(b.const, -tau1)
  # constraints #10
  b.const = c(b.const, tau2)
  
  
  # Build the LHS matrix
  ncol.A = d.x + d.x*n.obs + n.obs + d.x + d.f 
  I.dx = diag(d.x)
  I.df = diag(d.f)
  
  # const #1: beta upper bounds
  A1 = rbind(I.dx, -I.dx)
  A1 = cbind(A1, matrix(0, nrow=nrow(A1), ncol=ncol.A-d.x))
  
  # const #2: gamma upper bounds
  A2 = rbind(I.df, -I.df)
  A2 = cbind(matrix(0, nrow=nrow(A2), ncol=(ncol.A - d.f)), A2 )
  
  # const #3: delta.tilde upper bounds
  A3 = rbind(I.dx, -I.dx)
  A3 = cbind(matrix(0, nrow=nrow(A3), ncol=(ncol.A - d.x - d.f)), A3, matrix(0, nrow=nrow(A3), ncol=d.f) ) 
  
  # const #4: 0 <= l.tilde <= dt.tilde  
  # LHS is imposed by lower bounds later
  # RHS is l.tilde - dt.tilde <= 0 
  A4.1 = I.dx %x% diag(n.obs) 
  A4.2 = -rep(1,n.obs) %x% I.dx 
  A4 = cbind(matrix(0, nrow=nrow(A4.1), ncol=d.x), A4.1, matrix(0, nrow=nrow(A4.1), ncol=n.obs), A4.2, matrix(0, nrow=nrow(A4.1), ncol=d.f) )
  
  # const #4.5: 0 <= l.tilde <= dt.tilde  
  # LHS inequality: -l.tilde <= 0
  A4.11 = -I.dx %x% diag(n.obs) 
  A4.5 = cbind(matrix(0, nrow=nrow(A4.11), ncol=d.x), A4.11, matrix(0, nrow=nrow(A4.11), ncol=n.obs), matrix(0, nrow=nrow(A4.11), ncol=d.x) , matrix(0, nrow=nrow(A4.11), ncol=d.f) )
  
  
  # const #5: Left inequality of f_t'gm
  A5.1 = diag(M+2*eta)
  A5.2 = -f
  A5 = cbind(matrix(0, nrow=nrow(A5.1), ncol=d.x), matrix(0, nrow=nrow(A5.1), ncol=n.obs*d.x), A5.1, matrix(0, nrow=nrow(A5.1), ncol=d.x), A5.2 )
  
  # const #6: Right inequality of f_t'gm
  A6.1 = diag(-M) 
  A6.2 = f
  A6 = cbind(matrix(0, nrow=nrow(A6.1), ncol=d.x), matrix(0, nrow=nrow(A6.1), ncol=n.obs*d.x), A6.1, matrix(0, nrow=nrow(A6.1), ncol=d.x), A6.2 )
  
  # const #7: Right inequality sum l_t <= d_t * sum (U.dt - L.dt)
  sum.bnd = sum(U.dt - L.dt)
  A7.1 = diag(n.obs) %x% t(rep(1,d.x))
  A7.2 = (-sum.bnd)*diag(n.obs) 
  A7 = cbind(matrix(0, nrow=nrow(A7.1), ncol=d.x), A7.1, A7.2, matrix(0, nrow=nrow(A7.1), ncol=d.x), matrix(0, nrow=nrow(A7.1), ncol=d.f) )
  
  # const #7.5: Left inequality [- sum l_t] <= 0
  A7.3 = diag(n.obs) %x% t(rep(-1,d.x))
  A7.5 = cbind(matrix(0, nrow=nrow(A7.3), ncol=d.x), A7.3, matrix(0, nrow=nrow(A7.3), ncol=n.obs), matrix(0, nrow=nrow(A7.3), ncol=d.x), matrix(0, nrow=nrow(A7.3), ncol=d.f) )
  
  
  # const #8: Right inequality of sum (l - delta)
  A8.1 = - diag(n.obs) %x% t(rep(1,d.x))
  A8.2 = sum.bnd*diag(n.obs)
  A8.3 = rep(1,n.obs)%x% t(rep(1,d.x))
  A8 = cbind(matrix(0, nrow=nrow(A8.1), ncol=d.x), A8.1, A8.2, A8.3, matrix(0, nrow=nrow(A8.1), ncol=d.f) )
  
  # const #8.5: Left inequality of sum (l - delta)
  A8.11 = diag(n.obs) %x% t(rep(1,d.x))
  A8.33 = - rep(1,n.obs)%x% t(rep(1,d.x))
  A8.5 = cbind(matrix(0, nrow=nrow(A8.11), ncol=d.x), A8.11, matrix(0, nrow=nrow(A8.11), ncol=n.obs), A8.33, matrix(0, nrow=nrow(A8.11), ncol=d.f) )
  
  # const #9: Left inequality of sum d_t
  A9 = matrix(-rep(1/n.obs,n.obs), 1, n.obs)
  A9 = cbind(matrix(0, nrow=1, ncol=d.x), matrix(0, nrow=1, ncol=d.x*n.obs), A9, matrix(0, nrow=1, ncol=d.x),  matrix(0, nrow=1, ncol=d.f) )
  
  # const #10: Right inequality of sum d_t
  A10 = matrix(rep(1/n.obs,n.obs), 1, n.obs)
  A10 = cbind(matrix(0, nrow=1, ncol=d.x), matrix(0, nrow=1, ncol=d.x*n.obs), A10, matrix(0, nrow=1, ncol=d.x),  matrix(0, nrow=1, ncol=d.f) )
  
  A.const=rbind(A1,A2,A3,A4,A4.5, A5,A6,A7,A7.5, A8, A8.5, A9, A10)
  
  return(list(A.const=A.const,b.const=b.const))
}

# Function Name: build_constraint_selection
#
# This function constructs the whole linear constraints: Matrix of the LHS and a vector of RHS will be generated
#
# parameters = {bt,    l,   d,   dt,   gm1,   gm2,   e}
#               d.x  d.x*T  T    d.x  d.f-p    p     p
build_constraint_selection = function(L.bt, U.bt, L.dt, U.dt, L.gm1, U.gm1, L.gm2, U.gm2, L.gm3, U.gm3, L.p, U.p, M, eta, d.x, d.f, n.obs, f, tau1, tau2, p){
  
  # Build the RHS vector of constatins
  
  # constraints #1
  b.const = c(U.bt, -L.bt)
  # constraints #2
  b.const = c(b.const, U.gm1, -L.gm1)
  # constraints #2.1
  b.const = c(b.const, rep(0, p))
  # constraints #2.2
  b.const = c(b.const, rep(0, p))
  # constraints #2.3
  b.const = c(b.const, U.gm3, -L.gm3)
  # constraints #3
  b.const = c(b.const, (U.dt - L.dt), rep(0,d.x))
  # constraints #4
  b.const = c(b.const, rep(0,d.x*n.obs))
  # constraints #4.5
  b.const = c(b.const, rep(0,d.x*n.obs))
  # constraints #5
  b.const = c(b.const, M+eta)
  # constraints #6
  b.const = c(b.const, rep(0,n.obs))
  # constraints #7
  b.const = c(b.const, rep(0,n.obs))
  # constraints #7.5
  b.const = c(b.const, rep(0,n.obs))
  # constraints #8
  b.const = c(b.const, rep(sum(U.dt - L.dt),n.obs))
  # constraints #8.5
  b.const = c(b.const, rep(0,n.obs))
  # constraints #9
  b.const = c(b.const, -tau1)
  # constraints #10
  b.const = c(b.const, tau2)
  # constraints #11
  b.const = c(b.const, -L.p)
  # constraints #12
  b.const = c(b.const, U.p)
  
  
  # Build the LHS matrix
  ncol.A = d.x + d.x*n.obs + n.obs + d.x + d.f + p
  I.dx = diag(d.x)
  I.df = diag(d.f)
  d.f.gm1 = d.f -  p - 1
  I.df.gm1 = diag(d.f.gm1)
  I.df.gm2 = I.p = diag(p)
  
  # const #1: beta upper bounds
  A1 = rbind(I.dx, -I.dx)
  A1 = cbind(A1, matrix(0, nrow=nrow(A1), ncol=ncol.A-d.x))
  
  # const #2: gamma1 bounds
  A2 = rbind(I.df.gm1, -I.df.gm1)
  A2 = cbind(matrix(0, nrow=nrow(A2), ncol=(d.x+d.x*n.obs+n.obs+d.x)), A2, matrix(0, nrow=nrow(A2), ncol=2*p+1) )
  
  #------------------------------------------------------------------------------------------------
  #
  # constraints 2.1 -- 2.3 are added for model selection
  #
  #------------------------------------------------------------------------------------------------
  # const #2.1: gm2 - U.gm2*e <= 0
  I.U.gm2 = I.df.gm2
  diag(I.U.gm2) = U.gm2
  A2.1 = cbind(I.df.gm2, 0, -I.U.gm2)  # Order: gm2, gm3 (for constant term -1), e 
  A2.1 = cbind(matrix(0, nrow=nrow(A2.1), ncol=(d.x+d.x*n.obs+n.obs+d.x+d.f.gm1)), A2.1)
  
  # const #2.2: - gm2 + L.gm2*e <= 0
  I.L.gm2 = I.df.gm2
  diag(I.L.gm2) = L.gm2
  A2.2 = cbind(-I.df.gm2, 0, I.L.gm2)
  A2.2 = cbind(matrix(0, nrow=nrow(A2.2), ncol=(d.x+d.x*n.obs+n.obs+d.x+d.f.gm1)), A2.2)
  
  # const #2.3: gamma3 bounds
  A2.3 = rbind(1, -1)
  A2.3 = cbind(matrix(0, nrow=nrow(A2.3), ncol=(d.x+d.x*n.obs+n.obs+d.x+d.f.gm1+p)), A2.3, matrix(0, nrow=nrow(A2.3), ncol=p) )
  
  #------------------------------------------------------------------------------------------------
  
  # const #3: delta.tilde upper bounds
  A3 = rbind(I.dx, -I.dx)
  A3 = cbind(matrix(0, nrow=nrow(A3), ncol=(ncol.A - d.x - d.f)), A3, matrix(0, nrow=nrow(A3), ncol=d.f) ) 
  
  # const #4: 0 <= l.tilde <= dt.tilde  
  # LHS is imposed by lower bounds later
  # RHS is l.tilde - dt.tilde <= 0 
  A4.1 = I.dx %x% diag(n.obs) 
  A4.2 = -rep(1,n.obs) %x% I.dx 
  A4 = cbind(matrix(0, nrow=nrow(A4.1), ncol=d.x), A4.1, matrix(0, nrow=nrow(A4.1), ncol=n.obs), A4.2, matrix(0, nrow=nrow(A4.1), ncol=d.f+p) )
  
  # const #4.5: 0 <= l.tilde <= dt.tilde  
  # LHS inequality: -l.tilde <= 0
  A4.11 = -I.dx %x% diag(n.obs) 
  A4.5 = cbind(matrix(0, nrow=nrow(A4.11), ncol=d.x), A4.11, matrix(0, nrow=nrow(A4.11), ncol=n.obs), matrix(0, nrow=nrow(A4.11), ncol=d.x) , 
               matrix(0, nrow=nrow(A4.11), ncol=(d.f+p) ) )
  
  
  # const #5: Left inequality of f_t'gm
  A5.1 = diag(M+2*eta)
  A5.2 = -f
  A5 = cbind(matrix(0, nrow=nrow(A5.1), ncol=d.x), matrix(0, nrow=nrow(A5.1), ncol=n.obs*d.x), A5.1, matrix(0, nrow=nrow(A5.1), ncol=d.x), A5.2, matrix(0, nrow=nrow(A5.1), ncol=p))
  
  # const #6: Right inequality of f_t'gm
  A6.1 = diag(-M) 
  A6.2 = f
  A6 = cbind(matrix(0, nrow=nrow(A6.1), ncol=d.x), matrix(0, nrow=nrow(A6.1), ncol=n.obs*d.x), A6.1, matrix(0, nrow=nrow(A6.1), ncol=d.x), A6.2, matrix(0, nrow=nrow(A5.1), ncol=p) )
  
  # const #7: Right inequality sum l_t <= d_t * sum (U.dt - L.dt)
  sum.bnd = sum(U.dt - L.dt)
  A7.1 = diag(n.obs) %x% t(rep(1,d.x))
  A7.2 = (-sum.bnd)*diag(n.obs) 
  A7 = cbind(matrix(0, nrow=nrow(A7.1), ncol=d.x), A7.1, A7.2, matrix(0, nrow=nrow(A7.1), ncol=d.x), matrix(0, nrow=nrow(A7.1), ncol=d.f+p) )
  
  # const #7.5: Left inequality [- sum l_t] <= 0
  A7.3 = diag(n.obs) %x% t(rep(-1,d.x))
  A7.5 = cbind(matrix(0, nrow=nrow(A7.3), ncol=d.x), A7.3, matrix(0, nrow=nrow(A7.3), ncol=n.obs), matrix(0, nrow=nrow(A7.3), ncol=d.x), matrix(0, nrow=nrow(A7.3), ncol=d.f+p) )
  
  
  # const #8: Right inequality of sum (l - delta)
  A8.1 = - diag(n.obs) %x% t(rep(1,d.x))
  A8.2 = sum.bnd*diag(n.obs)
  A8.3 = rep(1,n.obs)%x% t(rep(1,d.x))
  A8 = cbind(matrix(0, nrow=nrow(A8.1), ncol=d.x), A8.1, A8.2, A8.3, matrix(0, nrow=nrow(A8.1), ncol=d.f+p) )
  
  # const #8.5: Left inequality of sum (l - delta)
  A8.11 = diag(n.obs) %x% t(rep(1,d.x))
  A8.33 = - rep(1,n.obs)%x% t(rep(1,d.x))
  A8.5 = cbind(matrix(0, nrow=nrow(A8.11), ncol=d.x), A8.11, matrix(0, nrow=nrow(A8.11), ncol=n.obs), A8.33, matrix(0, nrow=nrow(A8.11), ncol=d.f+p) )
  
  # const #9: Left inequality of sum d_t
  A9 = matrix(-rep(1/n.obs,n.obs), 1, n.obs)
  A9 = cbind(matrix(0, nrow=1, ncol=d.x), matrix(0, nrow=1, ncol=d.x*n.obs), A9, matrix(0, nrow=1, ncol=d.x),  matrix(0, nrow=1, ncol=d.f+p) )
  
  # const #10: Right inequality of sum d_t
  A10 = matrix(rep(1/n.obs,n.obs), 1, n.obs)
  A10 = cbind(matrix(0, nrow=1, ncol=d.x), matrix(0, nrow=1, ncol=d.x*n.obs), A10, matrix(0, nrow=1, ncol=d.x),  matrix(0, nrow=1, ncol=d.f+p) )
  
  #------------------------------------------------------------------------------------------------
  #
  # constraints 11 & 12 are added for model selection
  #
  #------------------------------------------------------------------------------------------------
  # const #11: Left inequality of sum e_t
  A11 = matrix(1, 1, p)
  A11 = cbind(matrix(0, nrow=nrow(A11), ncol=(ncol.A-p)), -A11)
  
  # const #12: Right inequality of sum e_t
  A12 = matrix(1, 1, p)
  A12 = cbind(matrix(0, nrow=nrow(A11), ncol=(ncol.A-p)), A12)
  #------------------------------------------------------------------------------------------------
  
  A.const=rbind(A1,A2,A2.1,A2.2,A2.3,A3,A4,A4.5, A5,A6,A7,A7.5, A8, A8.5, A9, A10, A11, A12)
  
  return(list(A.const=A.const,b.const=b.const))
}


# Function Name: estimate
#
# This function gives the estimation results
#
estimate <- function (y, x, f, Q.obj, L.obj, objcon, A.const, b.const, L.bt, L.gm, params) {
  
  # Call Library
  library("gurobi")
  
  model <- list()
  n.obs = length(y)
  d.x = ncol(x)
  d.f = ncol(f)
  
  # Quadratic objective function
  model$Q       = Q.obj
  model$obj     = L.obj
  model$objcon  = objcon
  
  # Linear constraints
  model$A       = A.const
  model$rhs     = b.const
  model$sense   = rep('<=', length(model$rhs))
  model$lb      = c(L.bt, rep(0,d.x*n.obs), rep(0,n.obs),  rep(0,d.x), L.gm)         #Done
  #model$lb      = c(rep(-10^10,d.x), rep(-10^10,d.x*n.obs), rep(-10^10,n.obs),  rep(-10^10,d.x), rep(-10^10,d.f))         #Done
  
  model$vtype   = c(rep('C', d.x), rep('C', d.x*n.obs), rep('B',n.obs), rep('C', d.x), rep('C', d.f))   
  model$modelsense = 'min'
  result <- gurobi(model, params=params)

  return(result)
}


# Function Name: estimate_selection
#
# This function gives the estimation results with model selection 
#
estimate_selection <- function (y, x, f, Q.obj, L.obj, objcon, A.const, b.const, L.bt, L.gm, params, p) {
  
  # Call Library
  library("gurobi")
  
  model <- list()
  n.obs = length(y)
  d.x = ncol(x)
  d.f = ncol(f)
  
  # Quadratic objective function
  model$Q       = Q.obj
  model$obj     = L.obj
  model$objcon  = objcon
  
  # Linear constraints
  model$A       = A.const
  model$rhs     = b.const
  model$sense   = rep('<=', length(model$rhs))
  model$lb      = c(L.bt, rep(0,d.x*n.obs), rep(0,n.obs),  rep(0,d.x), L.gm1, L.gm2, L.gm3, rep(0,p))
  
  model$vtype   = c(rep('C', d.x), rep('C', d.x*n.obs), rep('B',n.obs), rep('C', d.x), rep('C', d.f-p), rep('C', p), rep('B', p))   
  model$modelsense = 'min'
  result <- gurobi(model, params=params)
  
  return(result)
}



# Function Name: step1_grid
#
# This function estimates gm by the grid search method
#
step1_grid <- function(y, x, f, grid, eta=1e-6){
  n.grid = nrow(grid)
  dim.x = ncol(x)
  result = matrix(NA, nrow=n.grid, ncol=(1+2*dim.x))   # Collect ap.hat and the objective function value
  for (i in (1:n.grid)){
    f.index = as.numeric(f %*% t(grid[i,]))
    x.reg = cbind(x, x*(f.index>eta*0.99))
    m = lm(y~x.reg-1)
    result[i,] = c(sum(m$resid^2),coef(m))
    if (is.na(sum(coef(m)))) {
      result[i,1]=Inf
    }
  }
  opt = which.min(result[,1])
  ap.hat = result[opt,-1]
  gm.hat = as.numeric(grid[opt,])
  
  return(list(ap.hat=ap.hat, gm.hat=gm.hat))
  
}

# Function Name: estimate_bt_dt
#
# This function estimate bt and dt by OLS given an estimate for gm
#
estimate_bt_dt <- function(y, x, f, gm, eta=1e-6){
  T = length(y)                             # Sample size
  dim.bt = dim.dt = ncol(x)                 # The number of columns in x, size of Beta and Delta
  
  indicator = (f %*% gm >= eta*0.99)                # Construct the indicator function with gm given
  x.ind = x * matrix(indicator, T, dim.bt)  # Construct regressors for Delta
  reg = cbind(x, x.ind)                     # Construct the whole regressors for Beta and Delta
  m1 = lm(y~-1+reg)                         # OLS and save the model as m1
  
  bt.hat = m1$coef[1:dim.bt]
  dt.hat = m1$coef[-(1:dim.bt)]
  
  return(list(bt.hat=bt.hat, dt.hat=dt.hat, all.result=m1))
}

# Function Name: estimate_gm
# 
# This function estimates gm by MIO given estimates for bt and dt
#
estimate_gm <- function(y, x, f, bt, dt, A, b, M,  tau1, tau2, params, eta, ...){
  
  # Call necessary libraries
  library('gurobi')   # Gurobi library for MIO
  
  # Data dictionary: declare parameter values.
  n.obs = length(y)         # Sample size
  dim.gm = ncol(f)         # The dimension of gamma
  gm.hat = rep(NA,dim.gm)  # Maximizer of the problem
  dim.b = length(b)
  # do we need?  d.x = ncol(x)         # The dimension of delta
  
  
  # Declare a model 
  model <- list()
  
  # Model objective function: model$obj
  ind.dt.x = x %*% dt
  obj.d = (ind.dt.x^2 - 2 * ( y- ( x %*% bt ) ) * ind.dt.x)/n.obs
  
  model$obj = c(obj.d, rep(0,dim.gm))
  model$objcon = mean((y- x%*%bt)^2)
  
  #-------------------------------------------------------------------------------------------
  # Model constraints: matrix, rhs constants, and direction of inequalities for the constrains
  #------------------------------------------------------------------------------------------

  #-------------------------------------------------------
  # Constraints 1: re f'gm
  
  # Constaint 1.1: '(d.t-1)(M.t + 2 eta) + eta <= f'.t * gamma' reformulated into '(M.t +2eta) * d.t - f'.t * gamma <= M.t + eta'
  A1.1 = cbind(diag(as.vector(M+ 2*eta)), -f)
  b1.1 = M +eta
  
  # Constaint 1.2: 'f'.t * gamma <= d.t * M.t' reformulated into '- M.t d.t + f'.t * gamma <= 0'
  A1.2 = cbind(diag(as.vector(-M)), f)
  b1.2 = rep(0,n.obs)
  
  #-------------------------------------------------------
  
  # Constaint 2: Bounds for gm
  # A * gamma <= b
  zero.mat.3 = matrix( 0, nrow=dim.b, ncol=n.obs )
  A.2 = cbind( zero.mat.3, A )
  b.2 = b
  
  #-------------------------------------------------------
  # Constraints 3: re sum d_t
  
  # Constraint 3.1: Left inequality of sum d_t,    -sum d_t <= -tau1
  Ave.vec = rep(1,n.obs)/n.obs
  A3.1 = c(-Ave.vec, rep(0,dim.gm))
  b3.1 = -tau1
  
  # Constraint #5: Right inequality of sum d_t,   sum d_t <= tau2
  A3.2 = c(Ave.vec, rep(0,dim.gm))
  b3.2 = tau2
  
  model$A          = rbind(A1.1, A1.2, A.2, A3.1, A3.2)
  model$rhs        = c(b1.1, b1.2, b.2 , b3.1, b3.2)
  model$sense      = c(rep('<=',nrow(model$A)))
  
  # Other model parameter setting
  model$vtype      = c(rep('B',n.obs),rep('C',dim.gm))
  model$modelsense = "min"
  model$lb         = c(rep(0,n.obs),rep(-10^5,dim.gm))
  
  result = gurobi(model, params)
  
  # Return the estimate for gamma 
  return(list(obj=result$objval, sol=result$x, d.t=result$x[1:n.obs], gm=result$x[-(1:n.obs)], result=result))  
  
}



# Function Name: estimate_gm_selection
# 
# This function estimates gm and select the model by MIO given estimates for bt and dt
#
estimate_gm_selection <- function(y, x, f1, f2, bt, dt, L.gm1, U.gm1, L.gm2, U.gm2, L.gm3, U.gm3, M,  tau1, tau2, params, eta, ld, p){
  
  # Call necessary libraries
  library('gurobi')   # Gurobi library for MIO
  
  # Data dictionary: declare parameter values.
  n.obs = length(y)         # Sample size
  dim.gm = length(L.gm1)+length(L.gm2)+length(L.gm3)   # The dimension of gamma
  dim.d = n.obs
  gm.hat = rep(NA,dim.gm)  # Maximizer of the problem
  #dim.b = length(b)
  d.f.gm1 = dim.gm -  p - 1
  I.df.gm1 = diag(d.f.gm1)
  I.df.gm2 = I.p = diag(p)
  
  
  
  # Declare a model 
  model <- list()
  
  # Model objective function: model$obj
  # Parameters: (d, gamma, e)   Tx1, d.fx1, px1, respectively
  
  ind.dt.x = x %*% dt
  obj.d = (ind.dt.x^2 - 2 * ( y- ( x %*% bt ) ) * ind.dt.x)/n.obs
  
  model$obj = c(obj.d, rep(0,dim.gm), rep(ld,p))
  model$objcon = mean((y- x%*%bt)^2)
  
  # Model constraints: matrix, rhs constants, and direction of inequalities for the constrains
  #-------------------------------------------------------
  
  # Constraints 1: re f'gm
  
  # Constaint 1.1: '(d.t-1)(M.t + 2 eta) + eta <= f'.t * gamma' reformulated into '(M.t +2eta) * d.t - f'.t * gamma <= M.t + eta'
  A1.1 = cbind(diag(as.vector(M+ 2*eta)), -f)
  A1.1 = cbind(A1.1, matrix(0,nrow(A1.1),p))
  b1.1 = M +eta
  
  # Constaint 1.2: 'f'.t * gamma <= d.t * M.t' reformulated into '- M.t d.t + f'.t * gamma <= 0'
  A1.2 = cbind(diag(as.vector(-M)), f)
  A1.2 = cbind(A1.2, matrix(0,nrow(A1.2),p))
  b1.2 = rep(0,n.obs)
  
  #-------------------------------------------------------
  # Constaint 2: Bounds on gamma
  
  # constraint 2.1: bounds on gm1
  A2.1 = rbind(I.df.gm1, -I.df.gm1)
  A2.1 = cbind(matrix(0, nrow=nrow(A2.1), ncol=dim.d), A2.1, matrix(0, nrow=nrow(A2.1), ncol=2*p+1) )
  b2.1 = c(U.gm1, -L.gm1)
  
  # const #2.2: gm2 - U.gm2*e <= 0
  I.U.gm2 = I.df.gm2
  diag(I.U.gm2) = U.gm2
  A2.2 = cbind(I.df.gm2, 0, -I.U.gm2)  # Order: gm2, gm3 (for constant term -1), e 
  A2.2 = cbind(matrix(0, nrow=nrow(A2.2), ncol=dim.d+d.f.gm1), A2.2)
  b2.2 = rep(0,p)
  
  # const #2.3: - gm2 + L.gm2*e <= 0
  I.L.gm2 = I.df.gm2
  diag(I.L.gm2) = L.gm2
  A2.3 = cbind(-I.df.gm2, 0, I.L.gm2)
  A2.3 = cbind(matrix(0, nrow=nrow(A2.3), ncol=dim.d+d.f.gm1), A2.3)
  b2.3 = rep(0,p)
  
  # const #2.4: gamma3 bounds
  A2.4 = rbind(1, -1)
  A2.4 = cbind(matrix(0, nrow=nrow(A2.4), ncol=dim.d+d.f.gm1+p), A2.4, matrix(0, nrow=nrow(A2.4), ncol=p) )
  b2.4 = c(U.gm3, -L.gm3)
  
  #-------------------------------------------------------
  # Constraints 3: re sum d_t
  
  # Constraint 3.1: Left inequality of sum d_t,    -sum d_t <= -tau1
  Ave.vec = rep(1,n.obs)/n.obs
  A3.1 = c(-Ave.vec, rep(0,dim.gm+p))
  b3.1 = -tau1
  
  # Constraint #5: Right inequality of sum d_t,   sum d_t <= tau2
  A3.2 = c(Ave.vec, rep(0,dim.gm+p))
  b3.2 = tau2
  
  #-------------------------------------------------------
  # Constraints 4: sum e_t
  
  # Constraint 4.1: Left inequality of sum e_t,   -sum e_m <= - L.p
  A4.0 = matrix(1, 1, p)
  A4.1 = cbind(matrix(0, nrow=nrow(A4.0), ncol=dim.d+dim.gm), -A4.0)
  b4.1 = -L.p
  
  # const #12: Right inequality of sum e_t
  A4.2 = cbind(matrix(0, nrow=nrow(A4.0), ncol=dim.d+dim.gm), A4.0)
  b4.2 = U.p
  #-------------------------------------------------------
  
  
  model$A          = rbind( A1.1, A1.2, A2.1, A2.2, A2.3, A2.4, A3.1, A3.2, A4.1, A4.2)
  model$rhs        = c( b1.1, b1.2, b2.1, b2.2, b2.3, b2.4, b3.1, b3.2, b4.1, b4.2)
  model$sense      = c(rep('<=',nrow(model$A)))
  
  # Other model parameter setting
  model$vtype      = c(rep('B',n.obs),rep('C',dim.gm), rep('B',p))
  model$modelsense = "min"
  model$lb         = c(rep(0,n.obs),rep(-10^5,dim.gm), rep(0,p))
  
  result = gurobi(model, params)
  
  
  # Return the estimate for gamma 
  return(list(obj=result$objval, sol=result$x, d.t=result$x[1:n.obs], gm=result$x[(n.obs+1):(n.obs+dim.gm)], e=result$x[(n.obs+dim.gm+1):(n.obs+dim.gm+p)], result=result))  
  
}










# Function Name: get_objval
#
# This function evaluates the objective function value given data and parameter values
#
get_objval = function(y,x,f,bt,dt,gm, eta=1e-6) {
  n.obs=length(y)
  index = f%*%gm >= eta*0.99
  y.hat = x%*%bt + (x%*%dt) * index
  u.hat = y - y.hat
  objval = (t(u.hat) %*% u.hat)/n.obs
  
  return(objval)
}


# This function estimates an impurse-response functions using the local projection method in Jorda (2005)
#
#   y_(t+h) = shock_{t-1}*alpha + x_{t-1}*beta + + eps for h=0,1,2,...,
#
#   CALL: 
#         c_estimate_gm.R 
#         c_estimate_bt_dt.R 
#         c_get_m.R
#         c_l2norm.R
#
#   INPUT:
#         y: T x 1 
#         x: T x dim(x), it must include the constant term
#         s: T x dim(s), shock variables
#         h: scalar, time horinzeon of the impulse-reponse function
#         print: logical, print the regression results for all h, defaul=FALSE
#   
#   OUTPUT:
#         irf: time.horizon x dim(s)
#
# 

irf_local_projection <- function (y, x, s, h, print=F) {
  
  T = nrow(as.matrix(x))
  reg = cbind(s,x)
  #print(colnames(reg))
  y = as.matrix(y)
  irf = matrix(NA,h+1,ncol(s))
  
  for (i.h in (0:h)){
    sub.dep.y = y[((i.h+1):T),]
    sub.reg = reg[(1:(T-i.h)),]
    m = lm(sub.dep.y~sub.reg-1)
    if (print){
      cat('\n','--------------','\n','h=',i.h,'\n')
    }
    irf[i.h+1,] = coef(m)[paste('sub.reg',colnames(s),sep='')]
  }
  
  return(list(irf=irf))
}


