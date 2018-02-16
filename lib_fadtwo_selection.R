# File Name: lib_FADTWO.R
# This library contains all functions necessary for "FActor Driven TWO-regime regresssion (FADTWO)"
#
# Latest Update: 
#   2018-01-27    Original Code
# 


fadtwo <- function(y,x,f1, f2, method="joint", L.bt=-1e06, U.bt=1e06, L.dt=-1e06, U.dt=1e06, L.gm1=-1e06, U.gm1=1e06, L.gm2=-1e06, U.gm2=1e06, L.gm3=-1e06, U.gm3=1e06, 
                   L.p, U.p, tau1=0.05, tau2=0.95, eta=1e-6, params=list(OutputFlag=1, FeasibilityTol=1e-9),
                   grid=NULL, ld) {

  # Model: y = x'bt + x'dt * 1{f'gm > 0} + eps
  #
  # Input:
  #   y: outcome variable
  #   x: covariates
  #   f1: included factors
  #   f2: factors to be selected
  #   L.bt / U.bt: Lower and upper bounds for bt. The dim should be equal to the dim of bt
  #   L.dt / U.dt: Lower and upper bounds for dt
  #   L.gm / U.gm: Lower and upper bounds for gm
  #   tau1: Lower bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   tau2: Upper bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   eta: effective zero
  #   params: parameters for gurobi engine
  #
  x = as.matrix(x)
  f1 = as.matrix(f1)
  f2 = as.matrix(f2)
  f = cbind(f1,f2, -1)
  n.obs = nrow(x)     # Number of observations
  d.x = ncol(x)       # Dimension of regressors
  d.f = ncol(f)       # Dimension of factors
  p = ncol(f2)

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
    L.obj = get_L(x=x,y=y,L.dt=L.dt, d.f=d.f, ld=ld, p=p)
    Q.obj = get_Q(x=x, L.dt=L.dt, d.f=d.f, p=p)
    objcon = mean(y^2)
    
    # 2. Build the constraints
    const = build_constraint(L.bt=L.bt, U.bt=U.bt, L.dt=L.dt, U.dt=U.dt, L.gm1=L.gm1, U.gm=U.gm1, L.gm2=L.gm2, U.gm2=U.gm2, L.gm3=L.gm3, U.gm3=U.gm3, L.p=L.p, U.p=U.p, M=M, eta=eta, 
                             d.x=d.x, d.f=d.f, n.obs=n.obs, f=f, tau1=tau1, tau2=tau2, p=p)
    A.const = const$A.const
    b.const = const$b.const
    
    # 3. Estimate the model
    result = estimate(y=y, x=x, f=f, Q.obj=Q.obj, L.obj=L.obj, objcon=objcon, A.const=A.const, b.const=b.const, L.bt=L.bt, L.gm=L.gm, params=params, p=p)
    opt.par = result$x
    bt.hat = opt.par[1:d.x]
    l.hat = opt.par[(d.x+1):(d.x*(n.obs+1))]
    d.hat = opt.par[(d.x*(n.obs+1) +1):(d.x*(n.obs+1) + n.obs)]
    dt.tilde = opt.par[((d.x*(n.obs+1))+ n.obs + 1):( (d.x*(n.obs+1))+ n.obs  + d.x )]
    # Note that dt.tilde = dt.hat - L.dt.  
    dt.hat = dt.tilde + L.dt  
    gm1.hat = opt.par[((d.x*(n.obs+1))+n.obs  + d.x + 1):((d.x*(n.obs+1))+n.obs  + d.x  + d.f - p)]
    gm2.hat = opt.par[((d.x*(n.obs+1))+n.obs  + d.x  + d.f - p + 1):((d.x*(n.obs+1))+n.obs  + d.x  + d.f)]
    gm.hat = c(gm1.hat, gm2.hat)
    e.hat = opt.par[((d.x*(n.obs+1))+n.obs  + d.x  + d.f + 1):((d.x*(n.obs+1))+n.obs  + d.x  + d.f + p)]
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
    cat('Gamma1.hat   =', gm1.hat, '\n')
    cat('Gammm2.hat   =', gm2.hat, '\n')
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
    b.gm = c(U.gm1, U.gm2, -L.gm1, -L.gm2)
    M = rep(NA,n.obs)
    for (i in (1:n.obs)){
      M[i] = get_m(f=f[i,],A=A.gm,b=b.gm)$m
    } 
    
    # Estimate gm.hat by MIO and Update ap.hat
    step2_1.out = estimate_gm(y=y, x=x, f=f, bt=ap.hat.step1[c(1:d.x)], dt=ap.hat.step1[-c(1:d.x)], A=A.gm, b=b.gm, 
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
    cat('Sample size =', n.obs, '\n')
    cat('Obj val     =', objval, '\n')
    cat('Beta.hat    =', bt.hat, '\n')
    cat('Delta.hat   =', dt.hat, '\n')
    cat('Gamma.hat   =', gm.hat, '\n')
    cat('-----------------------------------------', '\n')
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
get_Q <-function(x, L.dt, d.f, p){
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
get_L <-function(x, y, L.dt, d.f, ld, p){
  
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
# parameters = {bt,    l,   d,   dt,   gm1,   gm2,   e}
#               d.x  d.x*T  T    d.x  d.f-p    p     p
build_constraint = function(L.bt, U.bt, L.dt, U.dt, L.gm1, U.gm1, L.gm2, U.gm2, L.gm3, U.gm3, L.p, U.p, M, eta, d.x, d.f, n.obs, f, tau1, tau2, p){
  
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
estimate <- function (y, x, f, Q.obj, L.obj, objcon, A.const, b.const, L.bt, L.gm, params, p) {
  
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
  
  # Model constraints: matrix, rhs constants, and direction of inequalities for the constrains
  # Constaint 1: '(d.t-1)(M.t + 2 eta) + eta <= f'.t * gamma' reformulated into '(M.t +2eta) * d.t - f'.t * gamma <= M.t + eta'
  Constraint.1 = cbind(diag(as.vector(M+ 2*eta)), -f)
  
  # Constaint 2: 'f'.t * gamma <= d.t * M.t' reformulated into '- M.t d.t + f'.t * gamma <= 0'
  Constraint.2 = cbind(diag(as.vector(-M)), f)
  
  # Constaint 3: A * gamma <= b
  zero.mat.3 = matrix( 0, nrow=dim.b, ncol=n.obs )
  Constraint.3 = cbind( zero.mat.3, A )
  
  # const #4: Left inequality of sum d_t
  A10.1 = rep(1,n.obs)/n.obs
  Constraint.4 = c(-A10.1, rep(0,dim.gm))
  
  # const #5: Right inequality of sum d_t
  Constraint.5 = c(-A10.1, rep(0,dim.gm))
  
  
  
  model$A          = rbind( Constraint.1, Constraint.2, Constraint.3, Constraint.4, Constraint.5)
  model$rhs        = c( M+eta, rep(0,n.obs), b , -tau1, tau2)
  model$sense      = c(rep('<=',nrow(model$A)))
  
  # Other model parameter setting
  model$vtype      = c(rep('B',n.obs),rep('C',dim.gm))
  model$modelsense = "min"
  model$lb         = c(rep(0,n.obs),rep(-10^5,dim.gm))
  
  result = gurobi(model, params)
  
  
  # Return the estimate for gamma 
  return(list(obj=result$objval, sol=result$x, d.t=result$x[1:n.obs], gm=result$x[-(1:n.obs)], result=result))  
  
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

