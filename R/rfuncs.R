#c1

#' Gradient-function for the SIR model
#' @param t Implicit argument for time
#' @param y  A vector with values for the states
#' @param parms A vector with parameter values for the SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 26, by=1/10)
#' paras  = c(mu = 0, N = 1, beta =  2, gamma = 1/2)
#' start = c(S=0.999, I=0.001, R = 0)
#' out=ode(y=start, times=times, func=sirmod, parms=paras)
#' @export
sirmod=function(t, y, parms){
   S=y[1]
   I=y[2]
   R=y[3]
   
   beta=parms["beta"]
   mu=parms["mu"]
   gamma=parms["gamma"]
   N=parms["N"]
   
   dS = mu * (N  - S)  - beta * S * I / N
   dI = beta * S * I / N - (mu + gamma) * I
   dR = gamma * I - mu * R
   res=c(dS, dI, dR)
   list(res)
 }

#c2

#' Negative log-likelihood function for the chain-binomial model
#' @param S0 a scalar with value for S0
#' @param beta a scalar with value for beta
#' @param I a vector incidence aggregated at serial interval
#' @return the negative log-likelhood for the model 
#' @examples
#' twoweek=rep(1:15, each=2)
#' niamey_cases1=sapply(split(niamey$cases_1[1:30], twoweek), sum)
#' llik.cb(S0=6500, beta=23, I=niamey_cases1)
#' @export
#' @importFrom stats dbinom
llik.cb = function(S0,beta,I){
    n = length(I)
    S = floor(S0-cumsum(I[-n]))
    p = 1-exp(-beta*(I[-n])/S0)  
    L = -sum(dbinom(I[-1],S,p,log=TRUE))  
    return(L)
}


#' Function to simulate the chain-binomial model
#' @param S0 a scalar with value for S0
#' @param beta a scalar with value for beta
#' @return A data-frame with time series of susceptibles and infecteds
#' @examples
#' sim=sim.cb(S0=6500, beta=23)
#' @export
#' @importFrom stats rbinom
sim.cb=function(S0, beta){
I=1
S=S0
i=1
while(!any(I==0)){
i=i+1
I[i]=rbinom(1, size=S[i-1], prob=1-exp(-beta*I[i-1]/S0))
S[i]=S[i-1]-I[i]
}
out=data.frame(S=S, I=I)
return(out)
}

#' Gradient-function for the chain-SIR model
#' @param t Implicit argument for time
#' @param logx  A vector with values for the log-states
#' @param params A vector with parameter values for the chain-SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 10, by=1/52)
#' paras2  = c(mu = 1/75, N = 1, beta =  625, gamma = 365/14, u=5)
#' xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001, paras2["u"]-1)), R = 0.0001))
#' out = as.data.frame(ode(xstart2, times, chainSIR, paras2))
#' @export
chainSIR=function(t, logx, params) {
    x = exp(logx)
    u = params["u"]
    S = x[1]
    I = x[2:(u + 1)]
    R = x[u + 2]
    
    with(as.list(params), {
        dS = mu * (N - S) - sum(beta * S * I)/N
        dI = rep(0, u)
        dI[1] = sum(beta * S * I)/N - (mu + u * gamma) * I[1]
        if (u > 1) {
            for (i in 2:u) {
                dI[i] = u * gamma * I[i - 1] - (mu + u * gamma) * 
                  I[i]
            }
        }
        dR = u * gamma * I[u] - mu * R
        res = c(dS/S, dI/I, dR/R)
        list(res)
    })
}

#c3

#' Auxillary function used by llik.pc 
#' @param a a vector with the ages 
#' @param up a vector with upper age-bracket cut-offs
#' @param foi a vector with FoI
#' @return A vector with FoIs matched to data
#' @seealso{llik.pc}
#' @export
integrandpc=function(a, up, foi){
wh=findInterval(a, sort(c(0,up)))
dur=diff(sort(c(0,up)))
inte=ifelse(wh==1, foi[1]*a, sum(foi[1:(wh-1)]*dur[1:(wh-1)])+foi[wh]*(a-up[wh-1]))
return(inte)
 }


#' Function to estimate parameters for the picewise-constant catalytic model
#'
#' This function uses binomial likelihoods to estimate the picewise-constant FoI model from age-incidence data
#'
#' @param par a vector with initial guesses 
#' @param age a vector with the ages 
#' @param num a vector with number infected by age
#' @param denom a vector with number tested by age
#' @param up a vector with upper age-bracket cut-offs
#' @return The negative log-likelihhod for a candidate piecewise constant catalytic model
#' @examples
#' x=c(1,4,8,12,18,24)
#' para=rep(.1,length(x))
#' \dontrun{optim(par=log(para),fn=loglikpc, age=rabbit$a, num=rabbit$inf, denom=rabbit$n, up=x)}
#' @export
#' @importFrom stats integrate
llik.pc = function(par, age, num, denom, up) {
ll = 0
for (i in 1:length(age)) {
p = 1 - exp(-integrandpc(a=age[i], up = up, foi = exp(par)))
ll = ll + dbinom(num[i], denom[i], p, log = T)
}
return(-ll)
}

#logspline not documented

#c4

#' Gradient-function for the SEIR model
#' @param t Implicit argument for time
#' @param y  A vector with values for the states
#' @param parms A vector with parameter values for the SEIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 10, by=1/120)
#' paras  = c(mu = 1/50, N = 1, beta =  1000, sigma = 365/8, gamma = 365/5)
#' start = c(S=0.06, E=0, I=0.001, R = 0.939)
#' out=ode(y=start, times=times, func=seirmod, parms=paras)
#' @export
seirmod=function(t, y, parms){
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]

  mu=parms["mu"]
  N=parms["N"]
  beta=parms["beta"]
  sigma=parms["sigma"]
  gamma=parms["gamma"]

  dS = mu * (N  - S)  - beta * S * I / N
  dE = beta * S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R 
  res=c(dS, dE, dI, dR)
  list(res)
}

#c4
#' Gradient-function for the forced SEIR model
#' @param t Implicit argument for time
#' @param y  A vector with values for the states
#' @param parms A vector with parameter values for the SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 10, by=1/120)
#' paras  = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2, sigma = 365/8, gamma = 365/5)
#' start = c(S=0.06, E=0, I=0.001, R = 0.939)
#' out=ode(y=start, times=times, func=seirmod2, parms=paras)
#' @export
seirmod2=function(t, y, parms){
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]

 with(as.list(parms),{
  dS = mu * (N  - S)  - beta0 * (1+beta1*cos(2*pi*t))* S * I / N
  dE = beta0 * (1+beta1*cos(2*pi * t))* S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res=c(dS, dE, dI, dR)
  list(res)
})
} 

#' Gradient-function for the age-structured SIR model with possibly heterogeneous mixin
#' @param t Implicit argument for time
#' @param logx  A vector with log-values for the log-states
#' @param parms A vector with parameter values for the age-structured SIR system
#' @return A list of gradients 
#' @examples
#' a=rep(1,4)
#' n=length(a)
#' betaM=matrix(1, ncol=4, nrow=4)
#' pars =list(N=1, gamma=365/14, mu=0.02, sigma=0.2, beta=500, betaM=betaM,p=rep(0,4), a=a)
#' xstart<-log(c(S=rep(0.099/n,n), I=rep(0.001/n,n), R=rep(0.9/n,n)))
#' times=seq(0,10,by=14/365)
#' out=as.data.frame(ode(xstart, times=times, func=siragemod, parms=pars))
#' @export
siragemod = function(t, logx,  parms){
  n=length(parms$a)
  xx = exp(logx)
        S = xx[1:n]
        I = xx[(n+1):(2*n)]
        R = xx[(2*n+1):(3*n)]
        
  with(as.list(parms), {
                phi = (beta*betaM%*%I)/N
                dS = c(mu,rep(0,n-1)) - (phi+a)*S + c(0,a[1:n-1]*S[1:n-1])*(1-p) - mu*S
                dI = phi*S + c(0,a[1:n-1]*I[1:n-1]) -(gamma+a)*I - mu*I
                dR =  c(0,a[1:n-1]*S[1:n-1])*p + c(0,a[1:n-1]*R[1:n-1]) + gamma*I - a*R - mu*R
                res = c(dS/S,dI/I,dR/R)
    list((res))
  })
}


#c6

#' Function to simulate the stocastic TSIR
#'
#' Function to simulate the stocastic TSIR assuming stochasticity in transmission and a Poisson birth-death process
#'
#' @param alpha the exponent on I
#' @param B the birth rate
#' @param beta the transmission rate
#' @param sdbeta the standard deviation on beta
#' @param S0 the initial susceptible fraction
#' @param I0 the initial number of infecteds
#' @param IT the length of simulation
#' @param N the population size
#' @return A list with time series of simulated infected and susceptible hosts
#' @examples
#' out = SimTsir()
#' @export
#' @importFrom stats rnorm
#' @importFrom stats rpois
SimTsir=function(alpha=0.97, B=2300, beta=25, sdbeta=0,
    S0 = 0.06, I0=180, IT=520, N=3.3E6){
    lambda = rep(NA, IT)
    I = rep(NA, IT)
    S = rep(NA, IT)
    I[1] = I0
    lambda[1] = I0
    S[1] = S0*N
    for(i in 2:IT) {
        lambda[i] = rnorm(1, mean=beta, sd=sdbeta) * I[i - 1]^alpha * S[i - 1] /N
        if(lambda[i]<0) {lambda[i]=0}
        I[i] = rpois(1, lambda[i])
        S[i] = S[i - 1] + B - I[i]
    }
    list(I = I, S = S)
}

#' Function to simulate the seasonally-forced TSIR
#'
#' Function to simulate the stocastic TSIR assuming stochasticity in transmission and a Poisson birth-death process
#'
#' @param beta the seasonal transmission coefficients
#' @param alpha the exponent on I
#' @param B a vector of Births (the length of which determines the length of the simulation)
#' @param N the population size
#' @param inits a list containing initial S and I
#' @param type an argument "det" or "stoc" that determines wheter a deterministic or stochastic simulation is done
#' @return A list with time series of simulated infected and susceptible hosts
#' @examples
#' \dontrun{see chapter 8 in book}
#' @export
SimTsir2=function(beta, alpha, B, N,  inits = list(Snull = 0, Inull = 0), type = "det"){
    type = charmatch(type, c("det", "stoc"), nomatch = NA)
    if(is.na(type))
        stop("method should be \"det\", \"stoc\"")
        
    IT = length(B)
    s = length(beta)
    lambda = rep(NA, IT)  
    I = rep(NA, IT)
    S = rep(NA, IT)
    
    I[1] = inits$Inull
    lambda[1] = inits$Inull
    S[1] = inits$Snull
    
    for(i in 2:IT) {
        lambda[i] = beta[((i - 2) %% s) + 1] * S[i - 1] * (I[i - 1]^alpha)/N
        if(type == 2) {
                I[i] = rpois(1, lambda[i])
            }
        if(type == 1) {
            I[i] = lambda[i]
        }
        S[i] =S[i - 1] + B[i] - I[i]
    }
    return(list(I = I, S = S))
}

#c7 

#' Gillespie exact algorithm
#'
#' Function simulating a dynamical system using the Gillespie exact algorithm
#'
#' @param rateqs a list with rate equations 
#' @param eventmatrix a matrix of changes in state variables associated with each event
#' @param parameters a vector of parameter values
#' @param initialvals a vector of initial values for the states
#' @param numevents number of events to be simulated
#' @return A data frame with simulated time series
#' @examples
#' rlist=c(quote(mu * (S+I+R)), quote(mu * S), quote(beta * S * I /(S+I+R)), 
#'  quote(mu * I), quote(gamma * I), quote(mu*R))
#' emat=matrix(c(1,0,0,-1,0,0,-1,1,0,0,-1,0,0,-1,1,0,0,-1),ncol=3, byrow=TRUE)
#' paras  = c(mu = 1, beta =  1000, gamma = 365/20)
#' inits = c(S=100, I=2, R=0)
#' sim=gillespie(rlist, emat, paras, inits, 100)
#' @export
#' @importFrom stats rexp
gillespie=function(rateqs, eventmatrix, parameters, initialvals, numevents){
res=data.frame(matrix(NA, ncol=length(initialvals)+1, nrow=numevents+1))
names(res)=c("time", names(initialvals))
res[1,]=c(0, initialvals)
for(i in 1:numevents){
rat=sapply(rateqs, eval, as.list(c(parameters, res[i,])))
res[i+1,1]=res[i,1]+rexp(1, sum(rat))
whichevent=sample(1:nrow(eventmatrix), 1, prob=rat)
res[i+1,-1]=res[i,-1]+eventmatrix[whichevent,]
}
return(res)
}


#' Gillespie tau-leap algorithm
#'
#' Function simulating a dynamical system using the Gillespie tau-leap approximation
#'
#' @param rateqs a list with rate equations 
#' @param eventmatrix a matrix of changes in state variables associated with each event
#' @param parameters a vector of parameter values
#' @param initialvals a vector of initial values for the states
#' @param deltaT the tau-leap time interval
#' @param endT the time length of simulation
#' @return A data frame with simulated time series
#' @examples
#' rlist2=c(quote(mu * (S+E+I+R)), quote(mu * S), quote(beta * S * I/(S+E+I+R)), 
#'  quote(mu*E), quote(sigma * E), quote(mu * I), quote(gamma * I), quote(mu*R))
#' emat2=matrix(c(1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1),
#' ncol=4, byrow=TRUE)
#' paras  = c(mu = 1, beta =  1000, sigma = 365/8, gamma = 365/5)
#' inits = c(S=999, E=0, I=1, R = 0)
#' sim2=tau(rlist2, emat2, paras, inits, 1/365, 1)
#' @export
tau=function(rateqs, eventmatrix, parameters, initialvals, deltaT, endT){
time=seq(0, endT, by=deltaT)
res=data.frame(matrix(NA, ncol=length(initialvals)+1, nrow=length(time)))
res[,1]=time
names(res)=c("time", names(initialvals))
res[1,]=c(0, initialvals)
for(i in 1:(length(time)-1)){
rat=sapply(rateqs, eval, as.list(c(parameters, res[i,])))
evts=rpois(1,  sum(rat)*deltaT)
if(evts>0){
whichevent=sample(1:nrow(eventmatrix), evts, prob=rat, replace=TRUE)
mt=rbind(eventmatrix[whichevent,], t(matrix(res[i,-1])))
#strange fix:
mt=matrix(as.numeric(mt), ncol=ncol(mt))
res[i+1,-1]=apply(mt,2,sum)
res[i+1, ][res[i+1,]<0]=0
}
else{
res[i+1,-1]=res[i,-1]
}}
return(res)
}

#' Function to predict efficacy of outbreak-response vaccination campaign
#' @param R reproductive ratio
#' @param day first day of ORV campaign
#' @param vaccine_efficacy Vaccine efficacy
#' @param target_vaccination fraction of population vaccinated during ORV campaign
#' @param intervention_length duration of ORV campaign
#' @param mtime length of simulation
#' @param LP length of latent period
#' @param IP length of infectious period
#' @param N initial susceptible population size
#' @return A list of gradients 
#' @examples
#' red1=retrospec(R=1.8, 161, vaccine_efficacy=0.85, target_vaccination=0.5, 
#'  intervention_length=10, mtime=250, LP=8, IP=5, N=16000)
#' 1-red1$redn
#' @export
retrospec<-function(R,day, vaccine_efficacy,target_vaccination,intervention_length, mtime, LP=7, IP=7, N=10000){
  steps<-1:mtime
        out<-matrix(NA,nrow=mtime, ncol=3)
  xstrt<-c(S=1-1/N,E=0,I=1/N,R=0,K=0)   #starting values
  beta<- R/IP         #transmission rate
  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
                P = 0, Dt = 0, T = Inf)
  outv<-as.data.frame(ode(xstrt,steps,sivmod,par))
  fsv<-max(outv$K)

  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
             P = target_vaccination, Dt = intervention_length, T = day)
outi<-as.data.frame(ode(xstrt,steps,sivmod,par))
  fsi<-max(outi$K)

         res<-list(redn=fsi/fsv)
  return(res)
}


#' Gradient-function for the SIR model with outbreak-response vaccination
#' @param t Implicit argument for time
#' @param x  A vector with values for the states
#' @param parms A vector with parameter values for the SIR system
#' @return A list of gradients 
#' @seealso \code{\link{retrospec}}
#' @export
sivmod<-function(t,x,parms){
    S<-x[1]
    E<-x[2]
    I<-x[3]
    R<-x[4]
    K<-x[5]
    with(as.list(parms),{
      Q<- ifelse(t<T | t>T+Dt,0,(-log(1-P)/Dt))
      dS<- -B*S*I-q*Q*S
      dE<- B*S*I-r*E
      dI<- r*E - g*I
      dR<- g*I+q*Q*S
      dK<-r*E
      res<-c(dS,dE,dI,dR,dK)
      list(res)
    })
  }



#c9

#' Gradient-function for Coyne et al's rabies model
#' @param t Implicit argument for time
#' @param logx A vector with values for the log-states
#' @param parms A vector with parameter values for the dynamical system
#' @return A list of gradients for the log system
#' @examples
#' require(deSolve)
#' times  = seq(0, 50, by=1/520)
#' paras  = c(gamma = 0.0397, b = 0.836, a = 1.34, sigma = 7.5, 
#' alpha = 66.36, beta = 33.25, c = 0, rho = 0.8)
#' start = log(c(X=12.69/2, H1=0.1, H2=0.1, Y = 0.1, I = 0.1))
#' out = as.data.frame(ode(start, times, coyne, paras))
#' @export
coyne=function(t, logx, parms){
  x=exp(logx)
  X=x[1]
  H1=x[2]
  H2=x[3]
  Y=x[4]
  I=x[5]
  N = sum(x)
  
  with(as.list(parms),{
  dX = a * (X + I) - beta * X * Y - gamma * N * X  - (b + c) * X
  dH1= rho * beta * X * Y  - gamma * N * H1  - (b + sigma + c) * H1
  dH2= (1-rho) * beta * X * Y  - gamma * N * H2  - (b + sigma + c) * H2
  dY = sigma * H1  - gamma * N * Y  - (b + alpha + c) * Y
  dI = sigma * H2 - gamma * N * I  - (b + c) * I
  res=c(dX/X, dH1/H1, dH2/H2, dY/Y, dI/I)
  list(res)
})
}

#c10

#' Function to do  Lyapunov exponent calculations from a TSIR simulation
#'
#' Function to do  Lyapunov exponent calculations from a TSIR simulation
#'
#' @param I a vector containg the time series of Is
#' @param S vector containg the time series of Ss
#' @param bt the seasonal transmission coefficients
#' @param alpha the exponent on I
#' @param N the population size
#' @return An object of class lyap with the lyapunov exponent, values for the Jacobians, parameters and data
#' @examples
#' \dontrun{see chapter 10 in book}
#' @export
TSIRlyap=function(I, S, alpha, bt, N){
  IT <- length(I)
  s <- length(bt)
  j11=rep(NA, IT)
  j12=rep(NA, IT)
  j21=rep(NA, IT)
  j22=rep(NA, IT)
  #initial unit vector
  J=matrix(c(1,0),ncol=1)
  #loop over the attractor
  for(i in 1:IT) {
  j11=1 -  bt[((i - 1) %% s) + 1] * I^alpha/N
  j12=-( bt[((i - 1) %% s) + 1] * S * (I^(alpha - 1) * alpha)/N)
  j21= bt[((i - 1) %% s) + 1] * I^alpha/N
  j22= bt[((i - 1) %% s) + 1] * S * (I^(alpha - 1) * alpha)/N
    J<-matrix(c(j11[i],j12[i],j21[i],j22[i]), ncol=2, byrow=TRUE)%*%J
      }
res=list(lyap=log(norm(J))/IT, j11=j11, j12=j12, j21=j21, j22=j22, I=I, S=S, alpha=alpha, bt=bt, N=N)
class(res)="lyap"
return(res)
}

#' Function to calculate the local Lyapunov exponents for the TSIR 
#'
#' Function to calculate the local Lyapunov exponents from an object of class \code{lyap}.
#'
#' @param x an object of class \code{lyap} (normally from a call to \code{TSIRlyap})
#' @param m number of forward iterations on the attractor
#' @return An object of class llyap with the local Lyapunov exponent and S-I data
#' @examples
#' \dontrun{see chapter 10 in book}
#' @export
TSIRllyap=function(x, m=1){
llyap=rep(NA, length(x$I))
for(i in 1:(length(x$I)-m)){
J=matrix(c(1,0,0,1), ncol=2)
for(k in 0:(m-1)){J = matrix(c(x$j11[(i+k)], x$j12[(i+k)], x$j21[(i+k)], x$j22[(i+k)]), ncol = 2, byrow=TRUE)%*%J}
llyap[i]=log(max(abs(eigen(J)$values)))/m
}
res=list(llyap=llyap, I=x$I, S=x$S)
class(res)="llyap"
return(res)
}

#' Gradient-function for the SIRWS model
#' @param t Implicit argument for time
#' @param logy  A vector with values for the log(states)
#' @param parms A vector with parameter values for the SIRWS system
#' @return A list of gradients (in log-coordinates)
#' @examples
#' require(deSolve)
#' times  = seq(0, 26, by=1/10)
#' paras  = c(mu = 1/70, p=0.2, N = 1, beta = 200, omega = 1/10, gamma = 17, kappa=30)
#' start = log(c(S=0.06, I=0.01, R=0.92, W = 0.01))
#' out = as.data.frame(ode(start, times, sirwmod, paras))
#' @export
sirwmod=function(t, logy, parms){
  y=exp(logy)
   S=y[1]
   I=y[2]
   R=y[3]
   W=y[4]

   with(as.list(parms),{
   dS = mu * (1-p) * N  - mu * S  - beta * S * I / N + 2*omega * W
   dI = beta * S * I / N - (mu + gamma) * I
   dR = gamma * I - mu * R - 2*omega * R +  kappa * beta * W * I / N + mu*p*N
   dW = 2*omega * R - kappa * beta * W * I / N - (2*omega +mu)* W 
   res=c(dS/S, dI/I, dR/R, dW/W)
   list(res)
 })
 }


#c11

#' Function to generate a ring lattice
#' @param N the number of nodes
#' @param K the number of neighbors to which each node is connected so degree = 2xK
#' @return An object of class CM (contact matrix)
#' @examples
#' cm=ringlattice(N=20,K=4)
#' @export
#' @importFrom stats toeplitz
ringlattice=function(N,K){
CM=toeplitz(c(0,rep(1,K),rep(0,N-2*K-1),rep(1,K)) )
    class(CM)="cm"
    return(CM)
}

#' Function to plot an object of class CM
#' @param x an object of class cm
#' @param ... other arguments 
#' @return A plot of the contract matrix
#' @examples
#' cm=ringlattice(N=20,K=4)
#' \dontrun{plot(cm)}
#' @export
#' @importFrom graphics symbols
#' @importFrom graphics segments
plot.cm=function(x, ...){
N=dim(x)[1]
theta=seq(0,2*pi,length=N+1)
x2=cos(theta[1:N])
y2=sin(theta[1:N])
symbols(x2,y2, fg=0, circles=rep(1, N), inches=0.1, bg=1, xlab="", ylab="")
segx1=as.vector(matrix(x2, ncol=length(x2), nrow=length(x2), byrow=TRUE))
segx2=as.vector(matrix(x2, ncol=length(x2), nrow=length(x2), byrow=FALSE))
segy1=as.vector(matrix(y2, ncol=length(x2), nrow=length(x2), byrow=TRUE))
segy2=as.vector(matrix(y2, ncol=length(x2), nrow=length(x2), byrow=FALSE))
segments(segx1,segy1, segx2, segy2, lty=as.vector(x))
}

#' Function to generate a Watts-Strogats network
#' @param N the number of nodes
#' @param K the number of neighbors to which each node is connected so degree = 2xK
#' @param Prw the rewiring probability
#' @return An object of class CM (contact matrix)
#' @examples
#' cm2=WattsStrogatz(N=20, K=4, Prw=.3)
#' @export
#' @importFrom stats runif
WattsStrogatz=function(N, K, Prw){
CM=ringlattice(N=N, K=K)
CMWS=CM
tri=CM[upper.tri(CM)]
Br=rbinom(length(tri),1,Prw)  # specify which edges to break
a=0
for(i in 1:(N-1)){                  
   for(j in (i+1):N){
  a=a+1               
  if(Br[a]==1 & CMWS[i,j]==1){ # if "break" is specified in Br matrix
    CMWS[i,j]=CMWS[j,i]=0 # break edge
    tmp=i
    tmp2=c(i, which(CMWS[i,]==1))             
    while(any(tmp2==tmp)){tmp=ceiling(N*runif(1))} # search new edge  
    CMWS[i,tmp]=CMWS[tmp,i]=1 # make new edge
    }
  }
}
class(CMWS)="cm"
return(CMWS)
}

#' Function to calculate the degree distribution for an object of class CM
#' @param object an object of class cm
#' @param plot if TRUE a barplot of the degree distribution is produced 
#' @param ... other arguments 
#' @return A plot of the contract matrix
#' @examples
#' cm=WattsStrogatz(N=20, K=4, Prw=.3)
#' summary(cm)
#' @export
#' @importFrom graphics barplot
summary.cm=function(object, plot=FALSE, ...){
  x=table(apply(object, 2, sum))
  res=data.frame(n=x)
  names(res)=c("degree", "freq")
  if(plot) barplot(x, xlab="degree")
    return(res)
}

#' Function to generate a Barabasi-Albert network
#' @param N the number of nodes
#' @param K the number of neighbors to which each node is connected so degree = 2xK
#' @return An object of class CM (contact matrix)
#' @examples
#' cm3=BarabasiAlbert(200, 4)
#' @export
BarabasiAlbert=function(N, K){
#https://en.wikipedia.org/wiki/Barabasi-Albert_model#Algorithm
CM=matrix(0, ncol=N, nrow=N)
CM[1,2]=1
CM[2,1]=1

for(i in 3:N){                  
  probs=apply(CM, 1, sum)
  link=unique(sample(c(1:N)[-i], size=min(c(K, i-1)), prob=probs[-i]))
    CM[i, link]=CM[link, i]=1
}
class(CM)="cm"
return(CM)
}

#' Function to simulate an epidemic on a network
#'
#' Function to simulate a stochastic (discrete time) Reed-Frost SIR model on a social network
#' 
#' @param CM a contact matrix
#' @param tau the transmission probabiliy
#' @param gamma the recovery probabiliy
#' @return An object of class netSIR with infectious status for each node through time
#' @examples
#' cm1=BarabasiAlbert(N=200,K=2)
#' sim1=NetworkSIR(cm1,.3,0.1)
#' summary(sim1)
#' \dontrun{plot(sim1)}
#' @aliases netSIR
#' @export
NetworkSIR=function(CM,tau,gamma){
#generate SIR epidemic on a network specified by the contact matrix 
#CM = contact matrix
#tau = probability of infection across an edge
#gamma = probability of removal per time step 
N=dim(CM)[1]
I=matrix(rep(0,N),nrow=N,ncol=1)   # initialize infecteds 
S=matrix(rep(1,N),nrow=N,ncol=1)  # initialize susceptibles 
R=matrix(rep(0,N),nrow=N,ncol=1)  # initialize removed 
I1=sample(1:N, size=1)
I[I1,1]=1 # initialize 1 random infected
S[I1,1]=0     

t=1
while(sum(I[,t-1])>0 | t==1){
    t=t+1
    infneigh=CM%*%I[,t-1]
    pinf=1-(1-tau)^infneigh
    newI=rbinom(N, S[,t-1], pinf)
    newR=rbinom(N, I[,t-1], gamma)

    nextS=S[,t-1]-newI
    nextI=I[,t-1]+newI-newR
    nextR=R[,t-1]+newR

    I=cbind(I, nextI)
    S=cbind(S, nextS)
    R=cbind(R, nextR)
    }
  
res=list(I=I,S=S,R=R)
class(res)="netSIR"
return(res)
}

#' Function to summarize a netSIR object
#' @param object an object of class netSIR
#' @return A data-frame with the time series of suceptible, infected and recovered individuals
#' @param ... other arguments 
#' @seealso
#' \code{\link{netSIR}}
#' @export
summary.netSIR=function(object, ...){
t=dim(object$S)[2]
S=apply(object$S,2,sum)
I=apply(object$I,2,sum)
R=apply(object$R,2,sum)
res=data.frame(S=S,I=I,R=R)
return(res)
}

#' Function to plot a netSIR object
#' @param x an object of class netSIR
#' @param ... other arguments 
#' @seealso
#' \code{\link{netSIR}}
#' @export
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom graphics lines
plot.netSIR=function(x, ...){
y=summary(x)  
plot(y$S, type="b", xlab="time", ylab="")
lines(y$I, type="b", col="red") 
lines(y$R, type="b", col="blue")  
    legend("right",
        legend=c("S", "I", "R"),
        lty=c(1,1, 1),
        pch=c(1,1, 1),
        col=c("black", "red", "blue"))
}

#' Function to calculate R0 from a contact matrix
#' @param CM an object of class CM
#' @param tau = probability of infection across an edge
#' @param gamma = probability of removal per time step 
#' @return the R0
#' @examples
#' cm1=BarabasiAlbert(N=200,K=2)
#' r0fun(cm1, 0.3, 0.1)
#' @export
r0fun=function(CM, tau, gamma){
x=apply(CM, 2, sum)
(tau/(tau+gamma))*(mean(x^2)-(mean(x)))/mean(x)
}

#c14

#' The Nicholson-Bailey model
#'
#' Function to simulate the Nicholson-Bailey Parasit-host model
#'
#' @param R the host repruductive rate
#' @param a the parasite search efficiency
#' @param T the length of simulation (number of time-steps)
#' @param H0 initial host numbers
#' @param P0 initial parasitoid numbers
#' @return A list of simulated Host and Parasitoid numbers
#' @examples
#' sim= NB(R=1.1,a=0.1)
#' @export
 NB = function(R, a, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   #we provide default parameters except for R and a
   H=rep(NA, T) #host series
   P=rep(NA, T) #parasitoid series

   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the parasitoid series

   for(t in 2:T){
     H[t] = R * H[t-1] * exp(- a * P[t-1])
     P[t] = R * H[t-1] * (1-exp(- a * P[t-1]))
   } #end of loop

   #the two vectors of results are stored in a "list"
   res= list(H = H, P = P) 
   return(res)
} 



#Shiny-apps
#' Launch a shiny-app simulating May's Parasitoid-host Model model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{May.app}
#' @export
#' @importFrom shiny renderPlot
May.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("May's Parasitoid-host Model"),
sidebarPanel(
sliderInput("R", "Growth rate (R):", 1.1,
              min = 1, max = 2, step=.01),
sliderInput("a", "Search efficiency (a):", 0.1,
              min = 0, max = .5),
sliderInput("k", "aggregation (k):", 1.5,
              min = 0.1, max = 3, step=0.1),
numericInput("P0", "Initial parasitoid:", 10,
              min = 1, max = 100),
numericInput("H0", "Initial host:", 20,
              min = 1, max = 100),
numericInput("Tmax", "Tmax:", 100,
              min = 1, max = 500)
),
mainPanel(tabsetPanel(
  tabPanel("Simulation", plotOutput("plot1", height = 500)),
     tabPanel("Details", 
    withMathJax(
         helpText("MODEL:"),
             helpText("Host $$H_t = R H_{t-1} (1 + a P_{t-1})^k$$"),
          helpText("Parasitoid $$P_t = R H_{t-1} (1-(1 + a P_{t-1})^k)$$"),
          helpText("REFERENCE: May RM (1978) Host-parasitoid systems in patchy 
            environments: a phenomenological model. J Anim Ecol 47: 833-843")
)
)
)
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
 NB = function(R, a, k, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   #we provide default parameters except for R and a

   H=rep(NA, T) #host series
   P=rep(NA, T) #parasitoid series

   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the host series

   for(t in 2:T){
     H[t] = R * H[t-1] * (1+ a * P[t-1])^(-k)
     P[t] = R * H[t-1] * (1-(1+ a * P[t-1])^(-k))
     if(P[t-1]==0) break
   } #end of loop

   #the two vectors of results are stored in a "list"
   res= data.frame(H = H, P = P)
 
   #the list is passed out of this function
   return(res)
} #end of function



  output$plot1 <- renderPlot({

    sim= NB(R=input$R, a=input$a, k=input$k, H0=input$H0, P0=input$P0, T=input$Tmax)
    time = 1:input$Tmax

    plot(time, sim$H, type= "b",xlab = "Generations", ylab = "Abundance", 
      ylim = range(sim, na.rm=TRUE))
    points(time, sim$P, type = "b", col = "red")
     legend("topleft",
        legend=c("H", "P"),
        lty=c(1,1),
        pch=c(1,1),
        col=c("black", "red"))
   })
  }
)

#' Launch a shiny-app simulating the seasonal SEIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{SEIR.app}
#' @export
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics curve
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics title
#' @importFrom deSolve ode
SEIR.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Seasonally forced SEIR"),
sidebarPanel(
sliderInput("beta0", "Transmission (yr^-1):", 1000,
              min = 0, max = 3000),
sliderInput("beta1", "Seasonality:", 0,
              min = 0, max = 1),
sliderInput("Ip", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("oneoversigma", "Latent period (days):", 8,
              min = 1, max = 100),
sliderInput("mu", "birth rate (per 1000):", 0.02,
              min = 0, max = .1),
sliderInput("T", "Time range:",
                  min = 0, max = 100, value = c(0,20)),
checkboxInput("lg", "un-Log", TRUE)
),
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
       tabPanel("Details", 
           withMathJax(
       helpText("MODEL:"),
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\frac{\\beta(t) I S}{N}$$"),
            helpText("Exposed $$\\frac{dE}{dt} = \\frac{\\beta(t) I S}{N} - (\\mu+\\sigma) E$$"),
            helpText("Infectious $$\\frac{dI}{dt} = \\sigma E - (\\mu+\\gamma) I$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\gamma I - \\mu R$$"),
           helpText("Seasonality $$\\beta(t) =  \\beta_0 (1 + \\beta_1 cos(2 \\pi t))$$"),
           helpText("Reproductive ratio $$R_0 =  \\frac{\\sigma}{\\sigma +\\mu} \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N}$$"),             
            helpText("REFERENCE: Earn DJD, Rohani P, Bolker BM, Grenfell BT (2000) A simple model for complex dynamical transitions in epidemics.
             Science 287: 667-670")
           ))
  
  )
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  seirmod2=function(t, x, params){
  S=x[1]
  E=x[2]
  I=x[3]
  R=x[4]

  mu=params["mu"]
  N=params["N"]
  beta0=params["beta0"]
  beta1=params["beta1"]
  sigma=params["sigma"]
  gamma=params["gamma"]

  dS = mu * (N  - S)  - beta0 * (1+beta1*cos(2*pi*t))* S * I / N
  dE = beta0 * (1+beta1*cos(2*pi * t))* S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res=c(dS, dE, dI, dR)
  list(res)
} 

#require(deSolve) 

  output$plot1 <- renderPlot({

  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = input$mu, N = 1, beta0 = input$beta0, beta1 = input$beta1, sigma = 365/input$oneoversigma, gamma = 365/input$Ip)
  xstart = c(S=0.06, E=0, I=0.001, R = 0.939)
  R0 = round(with(as.list(paras), sigma/(sigma+mu)*beta0/(gamma+mu)), 1)
out=ode(y=xstart,
  times=times,
  func=seirmod2,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

par(mar = c(5,5,2,5))
#lg=ifelse(input$lg==TRUE, "y", "")
plot(x=out$time[sel], y=out$I[sel], ylab="fraction", xlab="time", type="l",
ylim=range(out[sel,-c(1,2, 5)]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="red")
 lines(x=out$time, y=out$E, col="blue")
title(paste("R0=", R0))
# lines(x=out$time, y=out$S, col="green")
par(new=T)
plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA, 
    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""))
axis(side = 4, col="green")
mtext(side = 4, line = 4, "S", col="green")
  legend("right",
        legend=c("I", "E", "S"),
        lty=c(1,1,1),
         col=c("red", "blue", "green"))
   })
  
output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = input$mu, N = 1, beta0 = input$beta0, beta1 = input$beta1, sigma = 365/input$oneoversigma, gamma = 365/input$Ip)
  xstart = c(S=0.06, E=0, I=0.001, R = 0.939)
  R0 = with(as.list(paras), sigma/(sigma+mu)*beta0/(gamma+mu))
 
  out=ode(y=xstart,
  times=times,
  func=seirmod2,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(out$S[sel], out$I[sel], log=ifelse(input$lg==TRUE, "xy", ""), type="l", xlab="fraction susceptible", ylab="fraction infected")
  abline(v=1/R0, col="green")
  f4n=function(x){paras["mu"]*(1-x)/(paras["beta0"]*x)}
  curve(f4n, min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("I-socline", "S-isocline"),
        lty=c(1,1),
         col=c("red", "green"))
 
  })

  }
)

#' Launch a shiny-app simulating the SEIRS model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{SEIRS.app}
#' @export
#' @importFrom stats D
SEIRS.app=shinyApp(
# This creates the User Interface (UI)
ui <- pageWithSidebar(
headerPanel("SEIRS periodicity"),
sidebarPanel(
sliderInput("beta", "Transmission (yr^-1):", 500,
              min = 0, max = 3000),
sliderInput("oneoveromega", "Immune duration (years):", 4,
              min = 0, max = 100),
sliderInput("Ip", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("oneoversigma", "Latent period (days):", 8,
              min = 1, max = 100),
sliderInput("oneovermu", "Life expectancy (years):", 10,
              min = 1, max = 100),
sliderInput("T", "Time range:",
                  min = 0, max = 100, value = c(0,20)),
checkboxInput("lg", "un-Log", TRUE)
),

mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
       tabPanel("Equations", 
           withMathJax(
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\frac{\\beta I S}{N} + \\omega R$$"),
            helpText("Exposed $$\\frac{dE}{dt} = \\frac{\\beta I S}{N} - (\\mu+\\sigma) E$$"),
            helpText("Infecitous $$\\frac{dI}{dt} = \\sigma E - (\\mu+\\gamma) I$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\gamma I - \\mu R - \\omega R$$"),
           helpText("Reproductive ratio $$R_0 =  \\frac{\\sigma}{\\sigma +\\mu} \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N}$$")
          ))
  
  )
)
),

# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {
seirsmod=function(t, x, params){
  S=x[1]
  E=x[2]
  I=x[3]
  R=x[4]

  mu=params["mu"]
  beta=params["beta"]
  omega=params["omega"]
  sigma=params["sigma"]
  gamma=params["gamma"]

  dS = mu * (1  - S)  - beta * S * I + omega * R
  dE = beta * S * I - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R - omega * R
  res=c(dS, dE, dI, dR)
  list(res)
} 

  output$plot1 <- renderPlot({

  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = 1/input$oneovermu, beta =  input$beta, sigma = 365/input$oneoversigma, gamma = 365/input$Ip, omega=1/input$oneoveromega)
  xstart = c(S=0.539, E=0, I=0.001, R = 0.46)

  R0=with(as.list(paras),{
    beta*sigma/((mu+sigma)*(mu+gamma))
  })


Sstar=1/R0
Istar=paras["mu"]*(1-Sstar)/(paras["beta"]*Sstar - (paras["omega"]*paras["gamma"])/(paras["mu"]+paras["omega"]))
Estar=(paras["mu"]+paras["gamma"])*Istar/paras["sigma"]
Rstar=paras["gamma"]*Istar/(paras["mu"]+paras["omega"])

star=as.list(c(S=Sstar, E=Estar, I=Istar, R=Rstar, paras))
names(star)[1:4]=c("S", "E", "I", "R")

fns=list(quote(mu * (1  - S)  - beta * S * I  + omega * R), quote(beta * S * I - (mu + sigma) * E), quote(sigma * E - (mu + gamma) * I), quote(gamma * I - mu * R - omega * R))

aa1=as.vector(sapply(fns, D, "S"))
aa2=as.vector(sapply(fns, D, "E"))
aa3=as.vector(sapply(fns, D, "I"))
aa4=as.vector(sapply(fns, D, "R"))

JJ=matrix(c(sapply(aa1, eval, star), sapply(aa2, eval, star),sapply(aa3, eval, star),sapply(aa4, eval, star)), ncol=4)

EE=eigen(JJ)$values
WW=which.max(Im(EE))
rp=2*pi/Im(EE[WW])


out=ode(y=xstart,
  times=times,
  func=seirsmod,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

par(mar = c(5,5,2,5))
plot(x=out$time[sel], y=out$I[sel], ylab="fraction", xlab="time", type="l",
ylim=range(out[sel,-c(1,2, 5)]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="red")
 lines(x=out$time, y=out$E, col="blue")
 title(paste("R0=", round(R0, 1), ", Period=", round(rp,2)))

par(new=T)
plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA, 
    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""))
axis(side = 4, col="green")
mtext(side = 4, line = 4, "S", col="green")
  legend("right",
        legend=c("I", "E", "S"),
        lty=c(1,1,1),
         col=c("red", "blue", "green"))
   })
  
output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = 1/input$oneovermu, beta =  input$beta, sigma = 365/input$oneoversigma, gamma = 365/input$Ip, omega=1/input$oneoveromega)
  xstart = c(S=0.539, E=0, I=0.001, R = 0.46)

  R0=with(as.list(paras),{
    beta*sigma/((mu+sigma)*(mu+gamma))
  })
 

Sstar=1/R0
Istar=paras["mu"]*(1-Sstar)/(paras["beta"]*Sstar - (paras["omega"]*paras["gamma"])/(paras["mu"]+paras["omega"]))
Estar=(paras["mu"]+paras["gamma"])*Istar/paras["sigma"]
Rstar=paras["gamma"]*Istar/(paras["mu"]+paras["omega"])

star=as.list(c(S=Sstar, E=Estar, I=Istar, R=Rstar, paras))
names(star)[1:4]=c("S", "E", "I", "R")

fns=list(quote(mu * (1  - S)  - beta * S * I  + omega * R), quote(beta * S * I - (mu + sigma) * E), quote(sigma * E - (mu + gamma) * I), quote(gamma * I - mu * R - omega * R))

aa1=as.vector(sapply(fns, D, "S"))
aa2=as.vector(sapply(fns, D, "E"))
aa3=as.vector(sapply(fns, D, "I"))
aa4=as.vector(sapply(fns, D, "R"))

JJ=matrix(c(sapply(aa1, eval, star), sapply(aa2, eval, star),sapply(aa3, eval, star),sapply(aa4, eval, star)), ncol=4)

EE=eigen(JJ)$values
WW=which.max(Im(EE))
rp=2*pi/Im(EE[WW])


  out=ode(y=xstart,
  times=times,
  func=seirsmod,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(out$S[sel], out$I[sel], log=ifelse(input$lg==TRUE, "xy", ""), type="l", xlab="fraction susceptible", ylab="fraction infected")
  title(paste("R0=", round(R0, 1), ", Period=", round(rp,2)))
  abline(v=1/R0, col="green")
fffn=function(x){paras["mu"]*(1-x)/(paras["beta"]*x - (paras["omega"]*paras["gamma"])/(paras["mu"]+paras["omega"]))}
  curve(fffn, min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("I-socline", "S-isocline"),
        lty=c(1,1),
         col=c("red", "green"))
 
  })

  }
)



#' Launch a shiny-app simulating the SIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{SIR.app}
#' @export
#' @importFrom phaseR flowField
SIR.app=shinyApp(
# This creates the User Interface (UI)
ui <- pageWithSidebar(
headerPanel("The SIR model"),
sidebarPanel(
sliderInput("beta", "Transmission (yr^-1):", 300,
              min = 0, max = 1000),
sliderInput("infper", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("mu", "birth rate (/year):", 5,
              min = 0, max = 100),
sliderInput("T", "Time range:",
                  min = 0, max = 1, value = c(0,1))
),
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2", height = 500)),
      tabPanel("Equations", 
           withMathJax(
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\frac{\\beta I S}{N}$$"),
            helpText("Infecitous $$\\frac{dI}{dt} = \\frac{\\beta I S}{N} - (\\mu+\\sigma) I$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\gamma I - \\mu R$$"),
           helpText("Reproductive ratio $$R_0 =  \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N}$$")             
           ))
  
)   
  )
 )
,


# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {
  sirmod=function(t, x, parms){
    S=x[1]
    I=x[2]
    R=x[3]

    beta=parms["beta"]
    mu=parms["mu"]
    gamma=parms["gamma"]
    N=parms["N"]

    dS = mu * (N  - S)  - beta * S * I / N
    dI = beta * S * I / N - (mu + gamma) * I
    dR = gamma * I - mu * R
    res=c(dS, dI, dR)
    list(res)
  }

  output$plot1 <- renderPlot({
  times  = seq(0, input$T[2], by=1/1000)
  parms  = c(mu = input$mu, N = 1, beta =  input$beta, gamma =
    365/input$infper)
  start = c(S=0.999, I=0.001, R = 0)
  R0 = round(with(as.list(parms), beta/(gamma+mu)), 1)

  AA=with(as.list(parms), 1/(mu*(R0-1)))
  GG=with(as.list(parms), 1/(mu+gamma))
  rp=round(2*pi*sqrt(AA*GG),2)

  out=ode(y=start,
  times=times,
  func=sirmod,
  parms=parms)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(x=out$time[sel], y=out$S[sel], ylab="fraction", xlab="time", type="l",
  ylim=range(out[sel,-c(1,4)]))
  title(paste("R0=", R0, "Period=", rp))
  lines(x=out$time[sel], y=out$I[sel], col="red")
  lines(x=out$time[sel], y=out$R[sel], col="green")
  legend("right",
        legend=c("S", "I", "R"),
        lty=c(1,1,1),
         col=c("black", "red", "green"))
   })

  output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/1000)
  parms  = c(mu = input$mu, N = 1, beta =  input$beta, gamma =
    365/input$infper)
  start = c(S=0.999, I=0.001, R = 0)
  R0 = with(as.list(parms), beta/(gamma+mu))

AA=with(as.list(parms), 1/(mu*(R0-1)))
  GG=with(as.list(parms), 1/(mu+gamma))
  rp=round(2*pi*sqrt(AA*GG),2)
 
simod=function(t, y, parameters){
   S=y[1]
   I=y[2]

   beta=parameters["beta"]
   mu=parameters["mu"]
   gamma=parameters["gamma"]
   N=parameters["N"]
   
   dS = mu * (N  - S)  - beta * S * I / N
   dI = beta * S * I / N - (mu + gamma) * I
   res=c(dS, dI)
   list(res)
 }

  out=ode(y=start[-3],
  times=times,
  func=simod,
  parms=parms)

  out=as.data.frame(out)

  plot(x=out$S, y=out$I, xlab="Fraction suceptible", ylab="Fraction infected", type="l")
  title(paste("R0=", round(R0, 1), "Period=", rp))
  
fld=flowField(simod, xlim=range(out$S), ylim=range(out$I), 
parameters=parms, system="two.dim", add=TRUE,
ylab="I", xlab="S")
#cli=nullclines(simod, x.lim=c(0,.4), y.lim=c(0,.01), 
#parameters=parms, system="two.dim", add=TRUE, points=201)


  abline(v=1/R0, col="green")
  ffn=function(x){ parms["mu"]*(1-x)/(parms["beta"]*x)}
  curve(ffn, min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("I-socline", "S-isocline"),
        lty=c(1,1),
         col=c("red", "green"))

   })
  }
)

#' Launch a shiny-app simulating TSIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{TSIR.app}
#' @export
#' @importFrom stats spectrum
#' @importFrom polspline lspec
TSIR.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Simulating with TSIR"),
sidebarPanel(
sliderInput("alpha", "alpha:", 0.97,
              min = 0.8, max = 1),
sliderInput("beta", "beta:", 25,
              min = 0, max = 100),
sliderInput("B", "Births (B):", 2300,
              min = 0, max = 5000),
sliderInput("S0", "Fraction S:", 0.06,
              min = 0, max = 1),
numericInput("sdbeta", "sdbeta:", 3,
              min = 0, max = 10),
numericInput("N", "Popsize:", 3000000,
              min = 0, max = 100),
numericInput("I0", "initial I:", 100,
              min = 1, max = 100),
numericInput("IT", "Iterations:", 520,
              min = 1, max = 1000)
),
mainPanel(
  tabsetPanel(
      tabPanel("Simulation", plotOutput("plot1")), 
      tabPanel("Transfer function", plotOutput("plot2")),
      tabPanel("Details", 
           withMathJax(
            helpText("MODEL:"),
            helpText("Susceptible $$S_{t+1} = S_t - I_{t+1} + B$$"),
            helpText("Expected Inf $$\\lambda_{t+1} = \\frac{\\beta_t I_t^\\alpha S_t}{N}$$"),
           helpText("Infected $$I_{t+1} \\sim \\mbox{Poisson}(\\lambda_{t+1})$$"),
           helpText("Transmission $$\\beta_t \\sim \\mbox{Norm}(\\beta, \\mbox{sd}\\beta^2)$$"),             
            helpText("The transfer function is $$T(\\omega) = (\\vec{I}-e^{-\\imath \\omega} \\vec{J})^{-1} \\cdot \\vec{A}$$")
          )))   
  )),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  SimTsir=function(alpha, B, beta, sdbeta,
    S0, I0, IT, N){
    lambda = rep(NA, IT)
    I = rep(NA, IT)
    S = rep(NA, IT)
    I[1] = I0
    lambda[1] = I0
    S[1] = S0*N
    for(i in 2:IT) {
        lambda[i] = rnorm(1, mean=beta, sd=sdbeta) * I[i - 1]^alpha * S[i - 1] /N
        if(lambda[i]<0) {lambda[i]=0}
        I[i] = rpois(1, lambda[i])
        S[i] = S[i - 1] + B - I[i]
    }
    list(I = I, S = S)
}

  output$plot1 <- renderPlot({
  out = SimTsir(alpha=input$alpha, B=input$B, beta=input$beta, sdbeta=input$sdbeta, 
  S0=input$S0, I0=input$I0, IT=input$IT, N=input$N)
par(mfrow=c(1,2))  #This puts two plots side by side each other
plot(out$I, ylab="infected", xlab="time", type="b")
plot(out$S, out$I, ylab="infected", xlab="susceptible", type="b")
})

  output$plot2 <- renderPlot({
  out = SimTsir(alpha=input$alpha, B=input$B, beta=input$beta, sdbeta=input$sdbeta, 
  S0=input$S0, I0=input$I0, IT=input$IT, N=input$N)
Seq=expression(S-beta*S*I^alpha/N+B)
Ieq=expression(beta*S*I^alpha/N)
j11=D(Seq, "S")
j12=D(Seq, "I")
j21=D(Ieq, "S")
j22=D(Ieq, "I")
jj=c(j11, j12, j21,j22)

a1=D(Seq, "beta")
a2=D(Ieq, "beta")
aa=c(a1, a2)

paras  = c(B = input$B, beta =  input$beta, alpha = input$alpha, N=input$N)
eqs=sapply(c(quote(B^(1-alpha)*N/beta), quote(B)), eval, as.list(paras))

J=matrix(sapply(jj, eval, as.list(c(paras, c(S=eqs[1], I=eqs[2])))), 2, byrow=TRUE)
evs=eigen(J)$values
rp=2*pi/atan2(Im(evs[1]), Re(evs[1]))

A=matrix(sapply(aa, eval, as.list(c(paras, c(S=eqs[1], I=eqs[2])))), 2, byrow=TRUE)
Id=matrix(c(1,0,0,1),ncol=2)
wseq=seq(0,pi,length=500)

Fr=vector("list",500)  #set up empty list of matrices
# loop to fill those matrices with fourier trasform for the 500 values of w
for(i in 1:500){ 
Fr[[i]]=matrix(solve(Id-exp(1i*wseq[i])*J)%*%A,ncol=1)  #solve gives inverse
}

PS=matrix(NA,ncol=2,nrow=500,dimnames=list(1:500, c("S","I"))) 
#power spectra from real and imaginary parts of Fourier transform
for(i in 1:500){
PS[i,]=sqrt(Re(Fr[[i]])^2+Im(Fr[[i]])^2)
}

sfit=spectrum(out$I[-c(1:104)])
sfit2=lspec(out$I[-c(1:104)])
plot(wseq, PS[,2], type="l", 
	xlab="frequency (in radians)", ylab="amplitude", xlim=c(0,0.6))
title("Simulated spectrum (periodogram and log-spline) \n and T-fn prediction")
lines(pi*sfit$freq/0.5, max(PS[,2])*sfit$spec/max(sfit$spec), col=2)
par(new=TRUE)
plot(sfit2, col=3, xlim=c(0,0.6), axes=FALSE)
legend("topright", c("T-function", "periodogram", "log-spline"), lty=c(1,1,1), col=c(1,2,3))
   })
  }
)

#' Launch a shiny-app to study outbreak-response vaccination campaigns
#' @details
#' Launch app for details
#' @examples
#' \dontrun{orv.app}
#' @export
orv.app=shinyApp(ui=navbarPage("ORV",
  tabPanel("Intervention day",
    sidebarLayout(
    sidebarPanel(
sliderInput("vaccine_target1", "target:", 0.7,
              min = 0, max = 1),
sliderInput("vaccine_efficacy1", "efficacy:", 0.9,
              min = 0, max = 1),
numericInput("intervention_length1", "duration:", 14,
              min = 1, max = 100),
numericInput("N1", "N:", 1E5,
              min = 1E2, max = 1E6),
numericInput("mtime1", "endtime:", 120,
              min = 10, max = 3*365),
    sliderInput("R1", "R", 
              min = 1, max = 20, value=4),
    sliderInput("IP1", "Infectious period (days)", 5,
              min = 1, max = 100),
    sliderInput("LP1", "Latent period (days):", 8,
              min = 1, max = 100)),
    mainPanel(plotOutput("plot1"))
  )),
   tabPanel("R sensitivty",
    sidebarLayout(
    sidebarPanel(
    sliderInput("R2", "R", 
              min = 1, max = 20, value=4),
numericInput("pm2", "+/-:", 0.5,
              min = 1, max = 10),
sliderInput("vaccine_target2", "target:", 0.7,
              min = 0, max = 1),
sliderInput("vaccine_efficacy2", "efficacy:", 0.9,
              min = 0, max = 1),
numericInput("intervention_length2", "duration:", 14,
              min = 1, max = 100),
numericInput("N2", "N:", 1E5,
              min = 1E2, max = 1E6),
numericInput("mtime2", "endtime:", 120,
              min = 10, max = 3*365),
    sliderInput("IP2", "Infectious period (days)", 5,
              min = 1, max = 100),
    sliderInput("LP2", "Latent period (days):", 8,
              min = 1, max = 100)),
    mainPanel(plotOutput("plot2"))
  )),
   tabPanel("Duration sensitivity",
    sidebarLayout(
    sidebarPanel(
numericInput("intervention_length3", "duration:", 14,
              min = 1, max = 100),
numericInput("pm3", "+/-:", 7,
              min = 1, max = 21),
sliderInput("vaccine_target3", "target:", 0.7,
              min = 0, max = 1),
sliderInput("vaccine_efficacy3", "efficacy:", 0.9,
              min = 0, max = 1),
    sliderInput("R3", "R", 
              min = 1, max = 20, value=4),
numericInput("N3", "N:", 1E5,
              min = 1E2, max = 1E6),
numericInput("mtime3", "endtime:", 120,
              min = 10, max = 3*365),
    sliderInput("IP3", "Infectious period (days)", 5,
              min = 1, max = 100),
    sliderInput("LP3", "Latent period (days):", 8,
              min = 1, max = 100)),
    mainPanel(plotOutput("plot3"))
  )),
   tabPanel("Cover sensitivty",
    sidebarLayout(
    sidebarPanel(
sliderInput("vaccine_target4", "target:", 0.7,
              min = 0, max = 1),
numericInput("pm4", "+/-:", 0.1,
              min = 0, max = .9),
sliderInput("vaccine_efficacy4", "efficacy:", 0.9,
              min = 0, max = 1),
numericInput("intervention_length4", "duration:", 14,
              min = 1, max = 100),
    sliderInput("R4", "R", 
              min = 1, max = 20, value=4),
numericInput("N4", "N:", 1E5,
              min = 1E2, max = 1E6),
numericInput("mtime4", "endtime:", 120,
              min = 10, max = 3*365),
    sliderInput("IP4", "Infectious period (days)", 5,
              min = 1, max = 100),
    sliderInput("LP4", "Latent period (days):", 8,
              min = 1, max = 100)),
    mainPanel(plotOutput("plot4"))
  )),
   tabPanel("Retrospective analysis",
    sidebarLayout(
    sidebarPanel(
numericInput("day5", "Start day:", 60,
              min = 10, max = 3*365),
sliderInput("vaccine_target5", "target:", 0.7,
              min = 0, max = 1),
sliderInput("vaccine_efficacy5", "efficacy:", 0.9,
              min = 0, max = 1),
numericInput("intervention_length5", "duration:", 14,
              min = 1, max = 100),
sliderInput("R5", "R", 
              min = 1, max = 20, value=4),
numericInput("N5", "N:", 1E5,
              min = 1E2, max = 1E6),
numericInput("mtime5", "endtime:", 120,
              min = 10, max = 3*365),
    sliderInput("IP5", "Infectious period (days)", 5,
              min = 1, max = 100),
    sliderInput("LP5", "Latent period (days):", 8,
              min = 1, max = 100)),
    mainPanel(plotOutput("plot5"))
  )),
  tabPanel("Summary")
),

server=function(input, output){
######################################################
#SEIR model 
######################################################
simod<-function(t,x,parms){
################
#Parameters
#B = transmission rate
#1/r = latent period
#1/g = infectious period
#q = vaccine efficacy
#
#P = target coverage
#Dt = length of vaccination campaign
#T = Day of campaign start
#USAGE:
# times<-1:100
# xstrt<-c(S=.999,E=0,I=.001,R=0,K=0)
# par<-c(B=.5, r=1/7, g = 1/7, q = .8, P = 0, Dt = 10, T = 80)
# out<-as.data.frame(lsoda(xstrt,times,simod,par))
# plot(out$time,out$I,type="l")
#
#
# par<-c(B=.5, r=1/7, g = 1/7, q = .8, P = .99, Dt = 10, T = 50)
# out<-as.data.frame(lsoda(xstrt,times,simod,par))
# lines(out$time,out$I,col="red")

    S<-x[1]
    E<-x[2]
    I<-x[3]
    R<-x[4]
    K<-x[5]
    #
    with(as.list(parms),{
      Q<- ifelse(t<T | t>T+Dt,0,(-log(1-P)/Dt))
      dS<- -B*S*I-q*Q*S
      dE<- B*S*I-r*E
      dI<- r*E - g*I
      dR<- g*I+q*Q*S
      dK<-r*E
      res<-c(dS,dE,dI,dR,dK)
      list(res)
    })
  }
######################################################
  
#####################################################
#% Intervention
#####################################################
p_red<-function(R,vaccine_efficacy,target_vaccination,intervention_length, mtime=120, LP=7, IP=7, N=10000, step=1){
  steps<-(0:mtime)[seq(1,mtime,by=step)]
  p_red<-rep(NA,length(steps))
  xstrt<-c(S=1-1/N,E=0,I=1/N,R=0,K=0)   #starting values
  beta<- R/IP         #transmission rate
  t<-1
  for(i in 1:length(steps)){
    par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
                        P = target_vaccination, Dt = intervention_length, T = steps[i])
    out<-as.data.frame(lsoda(xstrt,steps,simod,par))
    p_red[t]<-out$K[dim(out)[1]]
    t<-t+1
    cat("step ", i,"of ",floor(mtime/step), ".\r")
  }
  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
                P = 0, Dt = 0, T = Inf)
  outv<-as.data.frame(lsoda(xstrt,steps,simod,par))
        #fs should really be the prediction with steps=Inf?
  fs<-max(out$K)
  res<-list(out=cbind(steps,p_red/max(p_red)),
                  R=R,
                  vaccine_efficacy=vaccine_efficacy,
                  target_vaccination=target_vaccination,
                  intervention_length=intervention_length,
                  mtime=mtime, LP=LP, IP=IP, N=N, step=step,
                  virgin=outv$I, vfs=fs)
  class(res)<-"p_red"
  return(res)
}
######################################################

#####################################################
#plotting pred objects
#####################################################
plot.p_red<-function(object){
      plot(object$out[,1],object$out[,2],type="l", xlab="First intervention day", ylab="% final epidemic", ylim=c(0,1))
      title(paste("target= ", round(100*object$target_vaccination,0), "% campaign = ", object$intervention_length,"d"))
par(new=TRUE)
      plot(object$out[,1], object$virgin, col='red', axes=FALSE, xlab="", ylab="", type="l")
      legend(x="topleft", legend=c("natural epidemic", "%final size"), col=c("red", "black"), lty=c(1,1))
}

#EX
#out<-p_red(R=4,vaccine_efficacy=.9,target_vaccination=.5,intervention_length=14, step=2)
#plot(out)

#####################################################
#% Sensitivity analysis on R
#####################################################
R_compare<-function(R=c(2,4,8),vaccine_efficacy=.9,target_vaccination, intervention_length,mtime=120, LP=7, IP=7, N=10000, step=7){
  out<-numeric(0)
  for(j in 1:length(R)){
    tmp<-p_red(R=R[j],vaccine_efficacy=vaccine_efficacy,
                target_vaccination,intervention_length, mtime=mtime, LP=LP, IP=IP, N=N, step=step)
    out<-cbind(out,tmp$out[,2])
  }
  res<-list(R=R, p_red=out, T=tmp$out[,1],
                 vaccine_efficacy=vaccine_efficacy,
                 target_vaccination=target_vaccination,
                 intervention_length=intervention_length,
                 mtime=mtime, LP=LP, IP=IP, N=N, step=step)
  class(res)<-"Rcomp"
  return(res)
}

#####################################################
#plotting pred objects
#####################################################
plot.Rcomp<-function(object){
        plot(NA,xlim=range(object$T),ylim=c(0,1), xlab="First intervention day", ylab="%, final epidemic")
        title(paste("% final size: target= ", round(100*object$target_vaccination,0), "% campaign = ", object$intervention_length,"d"))
  for(j in 1:length(object$R)){
    lines(object$T,object$p_red[,j],lty=j)
  }
      legend(x="right", legend=c(object$R), lty=1:length(object$R), title="R=")

}

#EX
#res<-R_compare(R=c(1.5, 2.5, 3.5), vaccine_efficacy=.9,target_vaccination=.5,intervention_length=7)
#plot(res)

#####################################################
#% Sensitivity analysis on length of intervention
#####################################################
Int_compare<-function(R,vaccine_efficacy,target_vaccination,intervention_length=c(7,10,14),mtime=120, LP=7, IP=7, N=10000, step=7){
  out<-numeric(0)
  for(j in 1:length(intervention_length)){
    tmp<-p_red(R=R,vaccine_efficacy=vaccine_efficacy,
                target_vaccination,intervention_length[j], mtime=mtime, LP=LP, IP=IP, N=N, step=step)
    out<-cbind(out,tmp$out[,2])
  }
  res<-list(R=R, p_red=out, T=tmp$out[,1],
                 vaccine_efficacy=vaccine_efficacy,
                 target_vaccination=target_vaccination,
                 intervention_length=intervention_length,
                 mtime=mtime, LP=LP, IP=IP, N=N, step=step)
  class(res)<-"Intcomp"
  return(res)
}

#####################################################
#plotting Intcomp objects
#####################################################
plot.Intcomp<-function(object){
        plot(NA,xlim=range(object$T),ylim=c(0,1), xlab="First intervention day", ylab="% final epidemic")
        title(paste("% final size: target= ", round(100*object$target_vaccination,0), "%, R = ", object$R))
  for(j in 1:length(object$intervention_length)){
    lines(object$T,object$p_red[,j],lty=j)
  }
      legend(x="right", legend=c(object$intervention_length), lty=1:length(object$intervention_length), title="Campaign D:")

}

#EX
#res<-Int_compare(intervention_length=c(7,10,14),R=3, vaccine_efficacy=.9, target_vaccination=.5)
#plot(res)

#####################################################
#% Sensitivity analysis on target Vaccination
#####################################################
Vacc_compare<-function(R,vaccine_efficacy,target_vaccination=c(.50,.70,.90),intervention_length=7,mtime=120, LP=7, IP=7, N=10000, step=7){
  out<-numeric(0)
  for(j in 1:length(target_vaccination)){
    tmp<-p_red(R=R,vaccine_efficacy=vaccine_efficacy,
                target_vaccination[j],intervention_length, mtime=mtime, LP=LP, IP=IP, N=N, step=step)
    out<-cbind(out,tmp$out[,2])
  }
  res<-list(R=R, p_red=out, T=tmp$out[,1],
                 vaccine_efficacy=vaccine_efficacy,
                 target_vaccination=target_vaccination,
                 intervention_length=intervention_length,
                 mtime=mtime, LP=LP, IP=IP, N=N, step=step)
  class(res)<-"Vacccomp"
  return(res)
}

#####################################################
#plotting Vaccomp objects
#####################################################
plot.Vacccomp<-function(object){
        plot(NA,xlim=range(object$T),ylim=c(0,1), xlab="First intervention day", ylab="% final epidemic")
        title(paste("% final size: R = ", object$R))
  for(j in 1:length(object$target_vaccination)){
    lines(object$T,object$p_red[,j],lty=j)
  }
      legend(x="right", legend=c(round(100*object$target_vaccination,0)), lty=1:length(object$target_vaccination), title="Target %")
}

#EX
#res<-Vacc_compare(target_vaccination=c(.50,.70,.90), R=4, vaccine_efficacy=.9,intervention_length=7)
#plot(res)

#####################################################
#Mortality analysis on R
#####################################################
M_red<-function(object, case_fatality){
        par(mfrow=c(2,1))
        plot(object)
  out<-(1-object$p_red)*object$N*case_fatality
  plot(NA,xlim=c(1,max(object$T)),ylim=c(0,max(out)), xlab="First intervention day", ylab="Extra survivors")
  for(j in 1:dim(out)[2]){
    lines(object$T, out[,j], lty=j)
  }
  title(paste("Reduced burden of mortality; N = ", object$N))
        par(mfrow=c(1,1))
}

#EX
#res<-Vacc_compare(target_vaccination=c(.50,.70,.90), R=4, vaccine_efficacy=.9,intervention_length=7)
#M_red(res, 0.5)


#####################################################
#Retrospective
#####################################################
retro<-function(R,day, vaccine_efficacy,target_vaccination,intervention_length, mtime=120, LP=7, IP=7, N=10000){
  steps<-1:mtime
        out<-matrix(NA,nrow=mtime, ncol=3)
  xstrt<-c(S=1-1/N,E=0,I=1/N,R=0,K=0)   #starting values
  beta<- R/IP         #transmission rate
  t<-1
  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
                P = 0, Dt = 0, T = Inf)
  outv<-as.data.frame(lsoda(xstrt,steps,simod,par))
        #fsv and fsi should really be with steps=Inf?
  fsv<-max(outv$K)

  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
             P = target_vaccination, Dt = intervention_length, T = day)

        outi<-as.data.frame(lsoda(xstrt,steps,simod,par))

  fsi<-max(outi$K)
        out[,1]<-steps
        out[,2]<-outv$I
        out[,3]<-outi$I
        res<-list(out=out,
                  redn=fsi/fsv,
                  R=R,
                  vaccine_efficacy=vaccine_efficacy,
                  target_vaccination=target_vaccination,
                  intervention_length=intervention_length,
                  mtime=mtime, LP=LP, IP=IP, N=N, day=day)
  class(res)<-"retro"
  return(res)
}
######################################################

#####################################################
#plotting pred objects
#####################################################
plot.retro<-function(object){
      plot(object$out[,1],object$out[,2],type="l", ylim=c(0,max(object$out[,2])), xlab='day', ylab='prevalence')
      polygon(c(object$day, object$day, object$day+object$intervention_length,
               object$day+object$intervention_length), c(-0.1,1,1,-.1), col="gray")
      lines(object$out[,1],object$out[,2])
      lines(object$out[,1], object$out[,3], col='red')
      title(paste("final size: ", round(100*(object$redn),1), "% (R=",
               object$R,", target=", 100*object$target_vaccination,"%)", sep=""))
      legend(x="topright", legend=c("natural epidemic", "w intervention"),
               col=c("black", "red"), lty=c(1,1))
      text(x=object$day+object$intervention_length, y=0, pos=4,
                labels=paste(object$intervention_length,
                "d intervention from", object$day))
}


output$plot1 <- renderPlot({
out<-p_red(R=input$R1,input$vaccine_efficacy1,input$vaccine_target1,input$intervention_length1, input$mtime1, input$LP1, input$IP1, input$N1, step=1)
plot(out)}
)

output$plot2 <- renderPlot({
R2=c(input$R2-input$pm2, input$R2, input$R2+input$pm2)
R2[R2<0]=0
out2=R_compare(R=R2, input$vaccine_efficacy2,input$vaccine_target2,input$intervention_length2, 
  input$mtime2, input$LP2, input$IP2, input$N2, step=1)
plot(out2)}
)

output$plot3 <- renderPlot({
il=c(input$intervention_length3-input$pm3,
input$intervention_length3, input$intervention_length3+input$pm3)
il[il<0]=0
out3=Int_compare(R=input$R3, input$vaccine_efficacy3,input$vaccine_target3,il,  
  input$mtime3, input$LP3, input$IP3, input$N3, step=1)
plot(out3)}
)

output$plot4 <- renderPlot({
vt=c(input$vaccine_target4-input$pm4,
input$vaccine_target4, input$vaccine_target4+input$pm4)
vt[vt<0]=0
vt[vt>1]=1
out4=Vacc_compare(R=input$R4, input$vaccine_efficacy4,vt,input$intervention_length4, 
 input$mtime4, input$LP4, input$IP4, input$N4, step=1)
plot(out4)}
)

output$plot5 <- renderPlot({
out5<-retro(R=input$R5, day=input$day5, input$vaccine_efficacy5,input$vaccine_target5,
  input$intervention_length5, input$mtime5, input$LP5, input$IP5, input$N5)
plot(out5)}
)

}
)
