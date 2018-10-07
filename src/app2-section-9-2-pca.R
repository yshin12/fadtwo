# Factor-driven two-regime regression
# An empirical exmple based on Hansen (1997)
#
# Last Update
#   2018-09-11 by Simon Lee  
#   2018-10-07
# Part of Bruce Hansen's R codes are used. 
# We thank him for posting his replication files on his web page.  

rm(list=ls())

library('ggplot2')
library('lmtest')
library('sandwich')
source('lib_fadtwo.R')

# ------------------------------------------------------------------------------------
# 
# Set the environment variables
#
#-------------------------------------------------------------------------------------


# Model selection
selection_method = 'no-selection' # 'l0' or 'no-selection'

# Estimation algorithm
method = 'joint' # 'joint' or 'iter'

# Number of additional factors
no_factors = 'pca' # 'pca', 'hansen', or 'all'

# Maximum iterations if method is 'iter'
K.bar = 2

# Tuning parameters for grid search
grid.type = 'fixed' # 'fixed' or 'random'
zeta = 0.5
grid.size = 10^6

# eta: the size of effective zero
eta = 1e-6

# Gurobi options: 
params <- list(OutputFlag=1, FeasibilityTol=1e-9, MIPGap=1e-4, TimeLimit=Inf)   
# OutputFlag: print out outputs
# FeasibilityTol: error allowance for the inequality constraints
# MIPGap: Stop if the gap is smaller than this bound. Default is 1e-4
# TimeLimit: Stop if the computation time reaches this bound. Default is infinity


# Minimum and maximum bounds for the propotion of each regime
tau1 = 0.15
tau2 = 0.85

# Number of lags in the regressor part
n.lags=12                                    

                    

# Specifiying the seed
set.seed(45462)
# ------------------------------------------------------------------------------------
# 
# End of setting the environment variables
#
#-------------------------------------------------------------------------------------


omit <- 0         # lags omitted from autoregression, if all included set omit=0       
name <- "UR"      # name of series                
p <- 12           # autoregressive order   


# Load Data into vector y #
u <- read.table("../data/LHMU-sample.dat")
c <- read.table("../data/LHMC-sample.dat")
factor1 <- read.table("../data/factor1-sample.dat")
factor2 <- read.table("../data/factor2-sample.dat")
u <- as.matrix(u)
c <- as.matrix(c)
factor1 <- as.matrix(factor1)
factor2 <- as.matrix(factor2)
s1 <- 3
n <- nrow(u)
y <- as.matrix((u[s1:n]/c[s1:n])*100)
factor1 <- factor1[s1:n]
factor2 <- factor2[s1:n]
n <- nrow(y)

# Taking First Differences #
dy <- as.matrix(y[2:n]-y[1:(n-1)])
factor1 <- factor1[1:(n-1)]
factor2 <- factor2[1:(n-1)]

# Create data matrix (lags, etc) #

n <- nrow(dy)

dat <- dy[(1+p):n]
factor1 <- factor1[(1+p):n]
factor2 <- factor2[(1+p):n]

names <- name
for (j in 1:p){
  if (sum(omit==j)==0){
    dat <- cbind(dat,dy[(1+p-j):(n-j)])
    if (j<10){ pn <- paste(c("0"),j,sep="")
    }else{ pn <- as.character(j)}
    namej <- paste(c("DY(t-"),pn,c(")"),sep="")
    names <- rbind(names,namej)    
  }  
}
xi <- seq(2,ncol(dat),by=1)
d <- 12
q <- y[(p+1):n] - y[(p+1-d):(n-d)]
dat = cbind(dat,q)
if (d<10) { pn <- paste(c("0"),d,sep="")
}else{ pn <- as.character(d)}
namej <- paste(c("Y*(t-"),pn,c(")"),sep="")
names <- rbind(names,namej)    



##############################################################

qi = ncol(dat)

hansen_factor = dat[,qi]
  
# dependent variable 
y = dat[,1]  
y = ts(y, start=c(1961,4), end=c(1996,7), frequency=12)   # 424 monthly observations from 1961.04 to 1886.07
dat.year <- time(y)  

# regressors

x = cbind(1,dat[,xi]) 


# Number of observations
n.obs = nrow(x)
# Number of covariates, dim(x)
d.x = ncol(x)
#-----------------------------------------------------------------------------------------------------------------------------
#
# Choose factors. 
#
#-----------------------------------------------------------------------------------------------------------------------------

if (no_factors == 'all'){
f = cbind(hansen_factor, factor1, -1)
}
if (no_factors == 'hansen'){
f = cbind(hansen_factor, -1)
}
if (no_factors == 'pca'){
f = cbind(factor1, -1)
}

# Number of factors, dim(f)
d.f = ncol(f)

# Decomposition of factors 
f1=as.matrix(f[,1])
if ((d.f > 2) && (selection_method == 'l0')){
f2=as.matrix(f[,2:(d.f-1)])
p = ncol(f2)  
}

# Run threshold regression with state = 1(y_{t-1} - y_{t-12}  >= 0.3020402)
    state.hansen = (hansen_factor > 0.3020402)
        x.hansen = cbind(x*as.numeric(1-state.hansen), x*as.numeric(state.hansen))
colnames(x.hansen) = c('_R1_const',paste('_R1_t',rep(1:12),sep='_'), '_R2_const',paste('_R2_t',rep(1:12),sep='_'))
        reg.hansen = lm(y~x.hansen-1)
 sigmahat.hansen = mean((y-fitted.values(reg.hansen))^2)
inference.hansen = coeftest(reg.hansen, vcov = vcovHC(reg.hansen, type = "HC3"))


output_text_file_name <- paste("../results/app2-section9-2",no_factors,".txt", sep="-")
sink(file = output_text_file_name, append = FALSE)
options(digits=3)
cat('----------------------------------------------------------------------------- \n')
cat("Estimation results using (L1.y-L.12) as the threshold variable (Hansen, 1997) \n")
cat("Coefficients are shown for two states (f1 < 0.3020402)) and (f1 > 0.3020402)), respectively \n")
print.table(round(inference.hansen,4))

sink()        

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
 
 
 if ((d.f > 2) && (selection_method == 'l0')){



L.p = 0
U.p = p
}
ld = sigmahat.hansen*log(n.obs)/n.obs     #BIC

time.start = proc.time()[3]



# Joint Estimation Algorithm
if (method == 'joint') {

if (selection_method == 'l0'){
    est.out=fadtwo_selection(y=y,x=x,f1=f1,f2=f2,method='joint',L.bt=L.bt,U.bt=U.bt,L.dt=L.dt,U.dt=U.dt,L.gm1=L.gm1,U.gm1=U.gm1,L.gm2=L.gm2,U.gm2=U.gm2,L.gm3=L.gm3,U.gm3=U.gm3,L.p=L.p, U.p=U.p, tau1=tau1,tau2=tau2, eta=eta,params=params,ld=ld)  
}
else {est.out=fadtwo(y=y,x=x,f=f,method='joint',L.bt=L.bt,U.bt=U.bt,L.dt=L.dt,U.dt=U.dt,L.gm=L.gm,U.gm=U.gm,tau1=tau1,tau2=tau2)
}    
}
  
# Iterative Estimation Algorithm
if (method == 'iter'){

  # Generate the grid points

if (grid.type == 'fixed'){   

  grid=gen_grid(option.grid='fixed', width=c(1,rep(zeta,d.f-1)), n.total=NULL, L.grid=c(L.gm1,rep(-Bnd.Const,d.f-1)), U.grid=c(U.gm1,rep(Bnd.Const,d.f-1)))
} else {
  grid=gen_grid(option.grid='random', width=NULL, n.total=grid.size, L.grid=c(L.gm1,rep(-Bnd.Const,d.f-1)), U.grid=c(U.gm1,rep(Bnd.Const,d.f-1)))
}  
  if (selection_method == 'l0'){
  est.out=fadtwo_selection(y=y,x=x,f1=f1,f2=f2,method='iter',L.gm1=L.gm1,U.gm1=U.gm1,L.gm2=L.gm2,U.gm2=U.gm2,L.gm3=L.gm3,U.gm3=U.gm3,L.p=L.p,U.p=U.p,tau1=tau1,tau2=tau2,eta=eta,params=params,grid=grid,max.iter=K.bar,ld=ld)
  }
  else {
  est.out=fadtwo(y=y,x=x,f=f,method='iter',L.gm=L.gm,U.gm=U.gm,tau1=tau1,tau2=tau2,grid=grid,max.iter=K.bar)
}
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

# Run threshold regression with the estimated state
    state.est = (f %*% gm.hat >= eta*.99)
        x.est = cbind(x*as.numeric(1-state.est), x*as.numeric(state.est))
colnames(x.est) = c('_R1_const',paste('_R1_t',rep(1:12),sep='_'), '_R2_const',paste('_R2_t',rep(1:12),sep='_'))
      reg.est = lm(y~x.est-1)
    resid.est = y-fitted.values(reg.est)
 sigmahat.est = mean(resid.est^2)
inference.est = coeftest(reg.est, vcov = vcovHC(reg.est, type = "HC3"))

sink(file = output_text_file_name, append = TRUE)
cat("------------------------------------------------------------------------------ \n")
cat("Estimation results with a vector of possible factors \n")
cat("Coefficients are shown for two states (f %*% gm.hat > 0) and (f %*% gm.hat < 0), respectively \n")
print.table(round(inference.est,4))

cat("sigmahat (Hansen) -- sigmahat (Est.) = ", sigmahat.hansen, ' -- ', sigmahat.est, '\n')


# Load NBER recession data #
        nber = read.table("../data/nber.dat")
        nber = ts(nber, start=c(1961,4), end=c(1996,7), frequency=12)   # 424 monthly observations from 1961.04 to 1886.07
  state.nber = (nber == 1)
match.hansen = 1-mean(abs(nber-as.numeric(state.hansen)))
   match.est = 1-mean(abs(nber-as.numeric(state.est)))


cat("Proportion of matches between NBER and Hansen = ", match.hansen, "\n")
cat("Proportion of matches between NBER and our estimates = ", match.est, "\n")

sink()



#-----------------------------------------------------------------------------------------------------------
#
# Graphical analysis
#
#----------------------------------------------------------------------------------------------------------


state.set = c('hansen','est','nber')

for (j_index in 1:3){
  if (state.set[j_index] == 'hansen'){state.y = state.hansen
  }
  if (state.set[j_index] == 'est'){state.y = state.est  
  }
  if (state.set[j_index] == 'nber'){state.y = state.nber  
  }
  
  # Draw different states on GNP
  u_figure_text_file_name <- paste("../results/app2-selection-9-2",state.set[j_index],no_factors, ".pdf", sep = "-")
  pdf(file=u_figure_text_file_name)
  recess = ggplot(data.frame(y, dat.year))+
    geom_line(mapping=aes_string(x="dat.year", y="y"))+
    annotate("rect", fill = "blue", alpha = 0.2, 
             xmin = dat.year[state.y], xmax = dat.year[state.y]+0.1,
             ymin = -Inf, ymax = Inf) +
    xlab('Year') +
    ylab('Unemployment Rates')
  print(recess)
  dev.off()
}


sink(file = output_text_file_name, append = TRUE)

time.end = proc.time()[3]
runtime = time.end - time.start
cat('Runtime     =', runtime, 'sec','\n')  
sink()  


save.image(file='app2-section-9-2-pca.RData')
