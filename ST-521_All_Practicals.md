---
title:"StochasticPracticals"
author:"VinaySanjayPatil"
date:"`rSys.Date()`"
output:
word_document:default
---

# Calculating limiting Distribution of MarkovChain
### 3.1ConsiderafiniteMarkovChain{Xn,n∈N}withstatespaceS={0,1,2,3}andTPMisgiven.FindthelimitingprobabilitiesofthisMarkovchain.

```{r}
library(expm)
P=matrix(c(2,1,1,0,
0,1,1,2,
1,2,1,0,
1,1,1,1),
ncol=4,byrow=TRUE)/4
#(a)Limitingdistributionvialinearequations
A=rbind(t(P)-diag(4),rep(1,4))
b=c(rep(0,4),1)
pi_eq=solve(t(A)%*%A,t(A)%*%b)
#(b)Powermethod
P50=P%^%50
P100=P%^%100
#(c)Differentinitialdistributions
alpha1=c(1,0,0,0)
alpha2=c(0.25,0.25,0.25,0.25)
alpha1%*%P100
alpha2%*%P100
pi_eq

```
### Conclusion:
#### 1.The chain is ergodic (irreducible + aperiodic), so it converges to a unique stationary distribution.

#### 2.The stationary distribution is: π=(0.2292, 0.3125, 0.25, 0.2083)

#### 3.No matter the starting state, after many transitions the probability of beingin each state stabilizes to these values.

#### 4.The fact that P^100 has identical rows equal to π confirms convergence.
#-----------------------------------------------------------------------------------------------------------------------------------------------

### 3.2 Consider a Markov chain {Xn,n∈N} with statespace S={0,1,2,3,4}. Suppose P04=1 and suppose that when the chain is in state i, i>0, the next state is equally likely to be any of the states 0,1,···,i−1. Find the limiting probabilities of this Markovchain.
```{r}
rm(list=ls())
library(expm)
P=matrix(c(0,0,0,0,1,
1,0,0,0,0,
1/2,1/2,0,0,0,
1/3,1/3,1/3,0,0,
1/4,1/4,1/4,1/4,0),ncol=5,byrow=TRUE)

#(b)Limitingdistribution
A=rbind(t(P)-diag(5),rep(1,5))
b=c(rep(0,5),1)
pi_eq=solve(t(A)%*%A,t(A)%*%b)

#(c)Powermethod
P100=P%^%100
#(d)DifferentInitialstates
alpha0=c(1,0,0,0,0)
alpha4=c(0,0,0,0,1)
alpha0%*%P100
alpha4%*%P100
pi_eq
```

# 4.Realization of Bienayme-Galton-Watson Branching Process

### 4.1 Let {Xn,n∈N} be a Bienayme–Galton–Watson branching process with offspring distribution P(ξ=0)=0.2, P(ξ=1)=0.3, P(ξ=2)=0.2, P(ξ=3)=0.3.
#### If X0=2, obtain realizations of X1,X2,...,X5.
#### (a) Generate one realization of the process.
#### (b) Repeat the realization using a different randomseed.
#### (c) Generate 10 independent realizations and observe the variability.
#### (d) Identify if extinction occurs in any realization.
```{r}
rm(list=ls())
set.seed(1)
z=c(0,1,2,3)
pz=c(0.2,0.3,0.2,0.3)
#(a)Singlerealization
x=numeric(6)
x[1]=2
for(i in 1:5)
x[i+1]=sum(sample(z,size=x[i],prob=pz,replace=TRUE))
x
#(b)Differentseed
set.seed(10)
x2=numeric(6)
x2[1]=2
for(i in 1:5){
x2[i+1]=sum(sample(z,size=x2[i],prob=pz,replace=TRUE))
return(x2[i+1])}
x2
#(c)Multiplerealizations
paths=matrix(0,nrow=10,ncol=6)
for(j in 1:10){
paths[j,1]=2
for(i in 1:5){
paths[j,i+1]=sum(sample(z,size=paths[j,i],prob=pz,replace=TRUE))
}
}
paths

```

### 4.2
```{r}
rm(list=ls())
set.seed(2)
z=c(0,1,2,3)
pz=c(0.25,0.35,0.3,0.1)
x=numeric(9)
x[1]=3
for(i in 1:8)
x[i+1]=sum(sample(z,size=x[i],prob=pz,replace=TRUE))
x
```

### 4.3
```{r}
set.seed(3)
z=c(0,1,2,3,4,5)
pz=c(0.1,0.2,0.15,0.2,0.25,0.1)
paths=matrix(0,nrow=10,ncol=8)
for(j in 1:10){
paths[j,1]=1

for(i in 1:7){
paths[j,i+1]=sum(sample(z,size=paths[j,i],prob=pz,replace=TRUE))
}
}
paths
```
# 5. Simulation of Poisson Process
### 5.1 Let {N(t),t≥0} be a Poisson process with rate λ=x.
### (a) Simulate the interarrival times and obtain the waiting times upto t=y.
### (b) Plot the sample path of N(t) for 0 ≤ t ≤ y.
### (c) Repeat the simulation using a different randomseed and compare the paths.
### Parameters:
### (i) x=2, y=5 (ii) x=4, y=9.5 (iii) x=9, y=7

```{r}
PoisProc=function(lambda,time){

inter=rexp(50,lambda)
arr=cumsum(inter)
arr=arr[arr<time]
n=length(arr)

t1=c(0,arr)
Nt=0:n

plot(t1,Nt,type="s",xlab="t",ylab="N(t)")
data.frame(t1,Nt)

}

set.seed(1)
PoisProc(2,5)

set.seed(2)
PoisProc(2,5)


```

### 5.2 Simulate the Poisson process at specified time points without simulating inter arrival times.
(a) Generate N(t) at the given time points.
(b) Repeat the simulation and compare outputs.
(c)Compare the simulated values with λt.
Parameters:
(a) t = 1.5, 2.2, 3.8, 7.5, 8.8, λ = 1
(b) t = 1.23, 2.21, 2.83, 6.05, 7.08, 17.8, λ = 1.5

```{r}
rm(list = ls())
SimuPoiss = function(lambda,t)
{
  n = length(t)
  Nt = numeric(n)
  Nt[1] = rpois(1, lambda*t[1])
  for ( i in 2:n){
    Nt[i] = Nt[i-1] + rpois(1, lambda*(t[i] - t[i-1]))
    
    data.frame(t, Nt, lambda_t = lambda*t)
  }

}

set.seed(1)
SimuPoiss(1,c(1.5,2.2,3.8,7.5,8.8))

```

### 5.3 Given that N(x) = y, simulate the Poisson process over 0 ≤ t ≤ x.
### (a) Generate a conditional realization of the process.
### (b) Plot the resulting path.
### (c) Repeat the experiment and compare paths.
### Parameters:
### (i) x = 6, y = 5 (ii) x = 7, y = 10 (iii) x = 9, y = 7

```{r}
rm(list = ls())
SimPoiss = function(n, t)
{
  arr = sort(runif(n,0,t))
  t1 = c (0 , arr )
  Nt = 0:n
  plot ( t1 , Nt , type="s" , xlab="t" ,     ylab="N( t ) " )
  data.frame(t1,Nt)
}
set.seed(1)
SimPoiss(5 ,6)
  
```


