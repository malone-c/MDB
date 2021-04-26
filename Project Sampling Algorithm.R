rm(list = ls())

require(mvtnorm)
require(timeR)

##### 1. LOAD IN THE DATA
murray = read.csv('data/murray_dummies.csv')

# Create the design matrix
intercept = rep(1, nrow(murray))
X = as.matrix(cbind(intercept, murray[, 3:ncol(murray)]))
y = log(murray$vol)

# Number of time periods
t = 2018 - 2005

# Number of panels
N = dim(murray)[1] / t

# Number of observations
n = t * N

# Fit a linear model to get starting values
OLS = lm(y ~ - 1 + as.matrix(X))

# Set the number of betas to estimate
m = ncol(X)

##### SET STARTING VALUES 
betas = rep(0, m)
sigma2 = summary(OLS)$sigma^2
rho = acf(OLS$residuals, plot = F, )$acf[2]

# Initialise the correlation matrix
DY.single = abs(outer((1:t), (1:t), "-"))
zero.matrix = matrix(0, nrow = N * t, ncol = N * t)
C.p = zero.matrix
corrs = rho^DY.single
for (i in seq(1, N * t, by = t)) {
	C.p[i : (i+t-1), i : (i + t-1)] = corrs
}

##### SET PRIOR HYPERPARAMETERS
Sigma.0 = diag(m)
beta.0 = rep(0, m)
nu.0 = 1
sigma2.0 = 1

# Get some useful numbers here to minimise redundant
# computation inside the algorithm.
a = (nu.0 + n) / 2
f.of.b = nu.0 * sigma2.0

inv.Sigma.0 = solve(Sigma.0)
inv.Sigma.0.times.beta.0 = inv.Sigma.0 %*% beta.0

##### SET S
S = 15000

#### SET DELTA FOR RHO
delta.rho = 0.09

### CREATE STORAGE OBJECT
store.beta = matrix(nrow = S, ncol = m)
store.rho = NULL
store.sigma2 = NULL

### CREATE ACCEPTANCE RATE OBJECTS
acc.rate.rho = 0


timer = createTimer()

{
timer$start("event1")

### START ALGORITHM
for (s in 1:S) {
	
	# Believe it or not, this is more than twice as fast as
	# just inverting C.p the usual way.
	inv.C.p = zero.matrix
	inv.corrs = solve(rho^DY.single)
	for (i in seq(1, N * t, by = t)) {
		inv.C.p[i : (i+t-1), i : (i + t-1)] = inv.corrs
	}
	
	
	# Gibbs for betas
	Sigma.n = solve((t(X) %*% inv.C.p %*% X) / sigma2 + inv.Sigma.0)
	b.n = Sigma.n %*% (((t(X) %*% inv.C.p %*% y) / sigma2) + inv.Sigma.0.times.beta.0)
	betas = t(rmvnorm(1, b.n, Sigma.n))
	
	# Gibbs for sigma2. Take y.minus.XB once to improve
	# efficiency.
	XB = X %*% betas
	y.minus.XB = y - XB
	SSR.rho = t(y.minus.XB) %*% inv.C.p %*% (y.minus.XB)
	sigma2 = 1 / rgamma(1, a, (f.of.b + SSR.rho) / 2)
	
	# Metropolis for rho
	rho.prop = abs(runif(1, rho - delta.rho, rho + delta.rho))
	rho.prop = min(rho.prop, 2 - rho.prop)
	
	# C.p.prop = zero.matrix
	# corrs.prop = rho.prop^DY.single
	# for (i in seq(1, N * t, by = t)) {
	# 	C.p.prop[i : (i+t-1), i : (i + t-1)] = corrs.prop
	# }
	
	# Do it with the covariance matrix instead to save time, so
	# we don't have to multiply C.p.prop (massive matrix) by
	# sigma2 when we compute densities in the next step.

	Sigma.prop = zero.matrix
	sigma.single.prop = sigma2 * (rho.prop^DY.single)
	for (i in seq(1, N * t, by = t)) {
		Sigma.prop[i : (i+t-1), i : (i + t-1)] = sigma.single.prop
	}
	
	r.rho.log = dmvnorm(y, XB, Sigma.prop, log = T) -
						dmvnorm(y, XB, sigma2 * C.p, log = T)
		
	
	# F. DRAW FROM U(0, 1)
			# (i) UPDATE PARAMETER
	    # (ii) UPDATE ACCEPTANCE RATE
	
	if (runif(1) < exp(r.rho.log)) {
		# Update rho and the correlation matrix
		rho = rho.prop

		# Update acceptance rate
		acc.rate.rho = acc.rate.rho + 1
	}
	

	
	store.beta[s, ] = betas
	store.sigma2 = c(store.sigma2, sigma2)
	store.rho = c(store.rho, rho)
	
	
		print(s)
	
}
timer$stop("event1")
}



# Diagnostics
# for (i in 1:m) {
# 	plot(store.beta[, i], type = 'l')
# }
# 
# plot(store.rho, type = 'l')
# acc.rate.rho / S
# 
# plot(store.sigma2, type = 'l')

# mdb = read.csv('data/murray.csv')
# lm = lm(log(vol) ~ price * crop + price * region + price * rainfall, mdb)
# length(lm$coefficients)
# step(lm)
# st = step(lm)
# summary(st)
# # Randomly remove half of the observations
# half = round(N / 2)
# sort(sample(1:N, half))


# plot(store.rho, type = 'l')
# plot(density(store.rho))
# 
# # Burn in period not required
# 
# 
# df = data.frame(beta = store.beta,
# 								rho = store.rho,
# 								sigma2 = store.sigma2,
# 								acc.rate.rho = acc.rate.rho)
# 
# index = seq(1, S, 5)
# coda::effectiveSize(store.rho)
# 
# plot(store.rho, type = 'l', xlim = c(0, 1000))
# abline(h = mean(store.rho), col = 'red')
# # Burn in period not required
# 
# 
# 
# plot(store.sigma2, type = 'l', xlim = c(0, 1000))
# abline(h = mean(store.sigma2), col = 'red')


coda::effectiveSize(store.beta)
head(murray)



