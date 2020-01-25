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
#' @return The negative log-likelihood for a candidate piecewise constant catalytic model
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

#' Gradient-function for the age-structured SIR model with possibly heterogeneous mixing
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

#' Function to simulate the stochastic TSIR
#'
#' Function to simulate the stochastic TSIR assuming stochasticity in transmission and a Poisson birth-death process
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
#' Function to simulate the stochastic TSIR assuming stochasticity in transmission and a Poisson birth-death process
#'
#' @param beta the seasonal transmission coefficients
#' @param alpha the exponent on I
#' @param B a vector of Births (the length of which determines the length of the simulation)
#' @param N the population size
#' @param inits a list containing initial S and I
#' @param type an argument "det" or "stoc" that determines whether a deterministic or stochastic simulation is done
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
#' @param I a vector containing the time series of Is
#' @param S vector containing the time series of Ss
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
#' @param K the number of neighbors to which each node is connected so degree = 2*K
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
#' @param plot if TRUE a bar plot of the degree distribution is produced 
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
#' @param K the number of neighbors to which each node is connected so degree = 2*K
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
#' @param tau the transmission probability
#' @param gamma the recovery probability
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
#' @return A data-frame with the time series of susceptible, infected and recovered individuals
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
#' Function to simulate the Nicholson-Bailey Parasitoid-host model
#'
#' @param R the host reproductive rate
#' @param a the parasitoid search efficiency
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

