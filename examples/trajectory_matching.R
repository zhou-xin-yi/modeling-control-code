# packages
require(ggplot2)
require(deSolve)
require(reshape)
require(bbmle)

## load data 
# Daily number of schoolboys confined to bed by influenza
baseURL <- "http://www.math.mcmaster.ca/bolker/eeid/data"
flu <- read.csv(url(paste(baseURL,"boarding_school_flu.csv",sep="/")))

ggplot(data=flu,mapping=aes(x=day,y=flu))+geom_point(size=2)+geom_line()

## ODE model
sibr.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  B <- x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  delta <- params["delta"]
  N <- 763
  ## now code the model equations
  dS.dt <- -beta*S*I/N
  dI.dt <- beta*S*I/N-gamma*I
  dB.dt <- gamma*I-delta*B
  ## combine results into a single vector
  dxdt <- c(dS.dt,dI.dt,dB.dt)
  ## return result as a list!
  list(dxdt)
}

params <- c(beta=2,gamma=1/3,delta=1/3)

times <- seq(from=0,to=15,by=1/4) ## returns a sequence
xstart <- c(S=762,I=1,B=0) ## initial conditions

out <- ode(
  func=sibr.model,
  y=xstart,
  times=times,
  parms=params
)

class(out)
head(out)
out <- as.data.frame(out)

# plot the result
out <- subset(out,select=c(I,B,time))
ggplot(data=melt(out,id.var="time"),
       mapping=aes(x=time,y=value,group=variable,color=variable))+
  geom_line()


## Incorporate the measurement error
times <- seq(from=0,to=14,by=1)

# negative binomial distribution
params <- c(beta=2,gamma=1/3,delta=1/3,N=763,p=0.3,k=1000)
out <- as.data.frame(
  ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=params
  )
)
within(
  out,
  C <- rnbinom(n=length(B),mu=params["p"]*B,size=params["k"])
) -> out
ggplot(data=out,mapping=aes(x=time,y=C))+geom_point()+geom_line()

## compute the likelihood of the data given the model and some parameters
# The following codes implement the negative log likelihood (nll) as a function of the model’s parameters
sibr.nll <- function (beta, gamma, delta, I.0, p, k) {
  times <- c(flu$day[1]-1,flu$day)
  ode.params <- c(beta=beta,gamma=gamma,delta=delta)
  xstart <- c(S=763-I.0,I=I.0,B=0)
  out <- ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=ode.params
  )
  ## 'out' is a matrix
  ll <- dnbinom(x=flu$flu,size=params["k"],mu=p*out[-1,"B"],log=TRUE)
  # return
  -sum(ll)
}

# The following codes plot the negative log likelihood as a function of beta
nll <- function (par) {
  sibr.nll(beta=par[1],gamma=1/3,delta=1/3,
           I.0=1,p=0.5,k=100)
}

betacurve <- data.frame(beta=seq(1/3,10,length=100)) # range of beta
within(betacurve,nll <- sapply(beta,nll)) -> betacurve
ggplot(data=betacurve,mapping=aes(x=beta,y=nll))+geom_line()

# find the exact value of the MLE for beta
fit <- optim(fn=nll,par=2,method="Brent",lower=1.5,upper=3)
fit$par # beta
fit$value # MLE


###
betacurve <- data.frame(beta=seq(1/3,10,length=100)) # range of beta
gamma <- vector("numeric", nrow(betacurve))
sibr.nll <- function (beta, gamma, delta, I.0, p, k) {
  times <- c(flu$day[1]-1,flu$day)
  ode.params <- c(beta=beta,gamma=gamma,delta=delta)
  xstart <- c(S=763-I.0,I=I.0,B=0)
  out <- ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=ode.params
  )
  ## 'out' is a matrix
  ll <- dnbinom(x=flu$flu,size=params["k"],mu=p*out[-1,"B"],log=TRUE)
  # return
  -sum(ll)
}
for (i in 1:length(betacurve)) {
  fit <- mle2(sibr.nll,parameters = c(betacurve$beta[i], ),method="Brent",lower = 0,upper = 1)
  gamma[i] <- coef(fit)[2]
}


ggplot(data=betacurve,mapping=aes(x=beta,y=nll))+geom_line()


# confidence interval
crit.lr <- pchisq(q=0.05,df=1,lower.tail=FALSE)
betacurve <- data.frame(beta=seq(2,3,length=100))
betacurve <- within(betacurve,nll <- sapply(beta,nll))
ggplot(data=betacurve,mapping=aes(x=beta,y=nll))+geom_line()+
  ylim(fit$value+c(0,10))+
  geom_vline(xintercept=fit$par,color='red')+
  geom_hline(yintercept=fit$value+crit.lr,color='blue')

# estimate both beta and gamma
nll <- function (par) {
  sibr.nll(beta=par[1],gamma=par[2],delta=1/3,
           I.0=1,p=0.5,k=100)
}
fit <- optim(fn=nll,par=c(2.6,1/3),method="Nelder-Mead")
fit$par # beta and gamma
fit$value # MLE

# mle2 command can isolate the MLE and compute profile likelihoods as well
fit <- mle2(sibr.nll,start=list(beta=2.5,gamma=0.7),
            method="Nelder-Mead",
            fixed=list(delta=1/3,k=100,p=0.5,I.0=1))
coef(fit)

pfit <- profile(fit)
plot(pfit)
confint(pfit)

# 3 parameters
fit <- mle2(sibr.nll,
            start=list(beta=1.5,gamma=1/3,p=0.5),
            method="L-BFGS-B",
            lower=c(0,0,0),
            upper=c(Inf,Inf,1),
            fixed=list(delta=1/3,k=100,I.0=1))
coef(fit)
pfit <- profile(fit)
#png("../plots/bbmle.png")
plot(pfit)
#dev.off()
confint(pfit)


# simulate the model
mle <- coef(fit)
times <- seq(from=0,to=14,by=1)
xstart <- c(S=763-mle["I.0"],I=mle["I.0"],B=0)
out <- as.data.frame(
  ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=mle
  )
)
within(
  subset(out,time>0),
  C <- rnbinom(n=length(B),mu=mle["p"]*B,size=mle["k"])
) -> out
ggplot(data=out,mapping=aes(x=time,y=C))+geom_line()+
  geom_point(data=flu,mapping=aes(x=day,y=flu))

