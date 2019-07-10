# Plotting t-distributions and standard normal curve
# Investigate CLT n < 30 if underlying distribution is not normal
curve(dt(x, df=5), xlim=c(4,-4), col="blue")
curve(dt(x, df=100), xlim=c(4,-4), col="red")
curve(dnorm(x, mean=0, sd=1), xlim=c(-4,4), col="green", add=T)

t_stat <- function (vector, mu) {
  t_statistic <- (mean(vector) - mu) / (sd(vector)/sqrt(length(vector)))
  if (is.infinite(t_statistic)) {
    return (NaN)
  } else {
    return (t_statistic)
  }
}

z_stat <- function (vector, mu, sd) {
  z_statistic <- (mean(vector) - mu) / (sd/sqrt(length(vector)))
  return (z_statistic)
}

# Sample s.d. is a biased estimator, except when n is large
n <- 10
x <-replicate(100000, rnorm(n,0,3))
s <- apply(x,2,sd)
mean(s)

# Simulation of z-test and t-test
n <- 5
x_small <-replicate(100000, rnorm(n,0,3))
small_t_distr <- apply(x_small,2,t_stat,0)
small_z_distr <- apply(x_small,2,z_stat,0,3)
range(small_z_distr)

hist(small_t_distr, breaks=seq(-30,30,0.2), xlim=c(-10,10), col="lightblue")
hist(small_z_distr, breaks=seq(-30,30,0.2), xlim=c(-10,10), col="pink")
par(new=T)
curve(dnorm, xlim=c(-10,10))
curve(dt(x,n-1), xlim=c(-10,10))

n <- 50
x <-replicate(100000, rnorm(n,0,3))
t_distr <- apply(x,2,t_stat,0)
z_distr <- apply(x,2,z_stat,0,3)

hist(t_distr, breaks=60, xlim=c(-10,10), col="lightblue")
hist(z_distr, breaks=60, xlim=c(-10,10), col="pink", add=T)
par(new=T)
curve(dnorm, xlim=c(-10,10))

# Z-test on exponential distribution
# Assumption of normality required for small n only!!!
# For large n, distribution need not be normal
n <- 5
lambda <- 3
x_exp <-replicate(10000, rexp(n, rate=lambda))

z_distr <- apply(x_exp, 2, z_stat, 1/lambda, sqrt(1/lambda^2))
hist(z_distr, breaks=50, xlim=c(-6,6))
par(new=T)
curve(dnorm, xlim=c(-6,6))

# T-test on exponential distribution
# Assumption of normality more strict!!!
t_distr <- apply(x_exp, 2, t_stat, 1/lambda)
hist(t_distr, breaks=100, xlim=c(-6,6))
par(new=T)
curve(dt(x,n-1), xlim=c(-6,6))

# CLT: Investigation of sampling distribution
n <- 3
lambda <- 3
# Exponential distribution with variance = 1/rate^2
# Investigate CLT Var(x_bar) = var(pop)/n
# True even for exponential!
x_exp <-replicate(10000, rexp(n, rate=lambda))
sample_mean <- apply(x_exp, 2, mean)
hist(sample_mean, breaks=30)
mean(sample_mean)
var(sample_mean)
clt_var <- 1/lambda/n

# CLT: Investigation of sampling distribution
# Investigating sampling distribution when underlying distribution is normal and n is small
# Sampling distribution still normal for small n
n <- 5
x_norm <-replicate(10000, rnorm(n))
sample_mean <- apply(x_norm, 2, mean)

hist(sample_mean, breaks=100, xlim=c(-2,2))
par(new=T)
curve(dnorm(x, 0, sqrt(1/n)), xlim=c(-2,2))

mean(sample_mean)
var(sample_mean)
1/n

x <- rnorm(10000)
A <- exp(x)
plot(density(log(A)),
     xlim = c(0,5))

# P-values of null distribution follow a uniform distribution
null_pvalue <- replicate(100000, t.test(rnorm(50), rnorm(50))$p.value)
hist(null_pvalue)

# Quantile-quantile plot
prob <- seq(0.001, 1 , 0.001)
x <- qnorm(prob)
y <- qnorm(prob, 0, 3)
y <- qt(prob, 10)
y <- sort(rt(1000, 3))

par(mfrow = c(1,1))
plot(x,y)
abline(a = 0, b = 1)

mix <- c(rnorm(500, mean = -3, sd = 1.5), rnorm(500, mean = 3, sd = 1.5))
hist(mix, breaks = 50)

length(x)
plot(x, sort(mix))
abline(a = 0, b = 1)
