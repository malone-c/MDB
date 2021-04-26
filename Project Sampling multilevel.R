rm(list = ls())

require(mvtnorm)
require(timeR)



##### 1. LOAD IN THE DATA
murray = read.csv('data/murray_dummies.csv')

# Design matrix
intercept = rep(1, nrow(murray))
X = as.matrix(cbind(intercept, murray[, 3:4]))
X[, 2] = log(murray[, 3])
y = log(murray$vol)

# Number of time periods
t = 2018 - 2005

# Number of panels
N = dim(murray)[1] / t

# Number of observations
n = t * N

# Number of betas to estimate
m = ncol(X)
p = 3

##### SET STARTING VALUES 
OLS = lm(y ~ - 1 + as.matrix(X))
sigma2 = summary(OLS)$sigma^2

# rho = acf(OLS$residuals, plot = F, )$acf[2]

# The ACF of the OLS residuals is ~0.84, but my model
# estimates the correlation coefficient to be closer to 0.7.
# So I set rho at 0.7, to hasten convergence.
rho = 0.7

Beta = matrix(0, nrow = N, ncol = p)
for (j in 1:N) {
	index = ((j-1) * t + 1) : (t * j)
	Beta[j, ] = lm(log(vol) ~ log(price) + rainfall, murray[index,])$coefficients
}

Sigma = cov(Beta)
inv.Sigma = solve(Sigma)
Theta = apply(Beta, 2, mean)

# Initialise the correlation matrix (and take inverse)
DY = abs(outer((1:t), (1:t), "-"))
corrs = rho^DY
inv.corrs = solve(corrs)

##### SET PRIOR HYPERPARAMETERS
mu.0 = Theta
Lambda.0 = cov(Beta)
inv.Lambda.0 = solve(Lambda.0)
nu.0 = 1
sigma2.0 = 1
S.0 = cov(Beta)
eta.0 = p + 2

# Get some useful numbers here to minimise redundant
# computation inside the algorithm.
a = (nu.0 + n) / 2
f.of.b = nu.0 * sigma2.0
zero.matrix = matrix(0, nrow = N * t, ncol = N * t)

##### SET S
S = 100000
thin = 10

#### SET DELTA FOR RHO
delta.rho = 0.05

### CREATE STORAGE OBJECTS
store.intercept = matrix(nrow = S / thin, ncol = N)
store.price = matrix(nrow = S / thin, ncol = N)
store.rainfall = matrix(nrow = S / thin, ncol = N)
store.theta = matrix(nrow = S / thin, ncol = p)
store.rho = NULL
store.sigma2 = NULL
acc.rate.rho = 0

# Set a timer to measure performance
timer = createTimer()

{
timer$start("event1")

### START ALGORITHM
for (s in 1:S) {
	
	# Gibbs for betas
	for (j in 1:N) {
		index = ((j-1) * t + 1) : (t * j)
		X.j = X[index, ]
		y.j = y[index]
		
		Sigma.n.j = solve((t(X.j) %*% inv.corrs %*% X.j) / sigma2 + inv.Sigma)
		
		mu.n.j = Sigma.n.j %*% (((t(X.j) %*% inv.corrs %*% y.j) / sigma2) + (inv.Sigma %*% (Theta)))
		Beta[j, ] = (rmvnorm(1, mu.n.j, Sigma.n.j))
	}
	
	# Gibbs for Theta (the means of the betas)
	Lambda.N = solve(inv.Lambda.0 + N * inv.Sigma)
	beta.mean = apply(Beta, 2, mean) 
	mu.N = Lambda.N %*% (inv.Lambda.0 %*% mu.0 + N * inv.Sigma %*% beta.mean)
	Theta = t(rmvnorm(1, mu.N, Lambda.N))

	# Gibbs for Sigma (the variance-covariance matrix of the
	# betas)
	Theta.matrix = matrix(Theta, N, p, byrow = T)
	S.theta = t(Beta - Theta.matrix) %*% (Beta - Theta.matrix)
	inv.Sigma = matrix(rWishart(1, eta.0 + N, solve(S.0 + S.theta)), 3, 3)
	
	# Gibbs for sigma2 (the variance of the y_i's)
	SSR.rho = 0
	for (j in 1:N) {
		index = ((j-1) * t + 1) : (t * j)

		y_j.minus.XB_j = y[index] - (X[index, ] %*% (Beta[j, ]))
		SSR.rho = SSR.rho + t(y_j.minus.XB_j) %*% inv.corrs %*% y_j.minus.XB_j
	}

	sigma2 = 1 / rgamma(1, a, (f.of.b + SSR.rho) / 2)
	
	# Get XB
	XB = NULL
	for (j in 1:N) {
		index = ((j-1) * t + 1) : (t * j)
		X.j.B.j = X[index, ] %*% Beta[j, ]
		XB = c(XB, X.j.B.j)
	}
	
	# Metropolis for rho (the temporal correlation coefficient).
	# Note that Sigma.prop is the variance-covariance matrix
	# for the proposed value of rho. It is not a proposed value
	# for Sigma, the variance-covariance of the betas, which is
	# handled above.
	rho.prop = abs(runif(1, rho - delta.rho, rho + delta.rho))
	rho.prop = min(rho.prop, 2 - rho.prop)
	
	C_p.prop = sigma2 * (rho.prop^DY)
	C_p.old = sigma2 * (corrs)
	
	numerator <- denominator <- 0
	for (j in 1:N) {
		index = ((j-1) * t + 1) : (t * j)
		
		numerator = numerator + dmvnorm(y[index], XB[index], C_p.prop, log = T)
		denominator = denominator + dmvnorm(y[index], XB[index], C_p.old, log = T)
	}
	
	r.rho.log = numerator - denominator
		
	if (runif(1) < exp(r.rho.log)) {
		rho = rho.prop
		acc.rate.rho = acc.rate.rho + 1
		
		# Perform matrix computation in here so that we only
		# compute new matrices (expensive) when we accept a new rho
		corrs = rho^DY
		inv.corrs = solve(corrs)
	}
	

	# Print the current iteration number to the console, for
	# monitoring purposes
	if (s %% thin == 1) {
		print(s)
		
		# Store all drawn parameters from this iteration
		# Each column corresponds to a group
		store.intercept[s / thin, ] = Beta[, 1]
		store.price[s / thin, ] = Beta[, 2]
		store.rainfall[s / thin, ] = Beta[, 3]
		store.theta[s / thin, ] = t(Theta)
		store.sigma2 = c(store.sigma2, sigma2)
		store.rho = c(store.rho, rho)
	}
}
	
# Stop the timer
timer$stop("event1")

# Hit the alarm so I know it's done
beepr::beep()
}



# Use these to figure out which crops and regions ahve
# highest and lowest elasticity
# crops.in.order = crop[seq(1, n - t, by = t)]
# region.in.order = region[seq(1, n - t, by = t)]
# 
# df = list(intercept = store.intercept,
# 								price = store.price,
# 								rainfall = store.rainfall,
# 								theta = store.theta,
# 								sigma2 = store.sigma2[1:10000],
# 					rho = store.rho[1:10000],
# 					acc.rate.rho = acc.rate.rho)
# write.csv(df, '100000draws_random_coef_FINAL.csv')

es = NULL
for (i in 1:50){
	es = c(es, coda::effectiveSize(store.rho[seq(1, S, by = i)]))
}
es
# sort(coda::effectiveSize(as.data.frame(df)))
# # 
# baba = read.csv('20000draws_random_coef.csv')

write.csv(store.price, 'FINAL_price_coefficients.csv')


# Diagnostics
for (i in 1:N) {
	plot(store.price[, i], type = 'l')
}

plot(store.rho, type = 'l')
acc.rate.rho / S

plot(store.sigma2, type = 'l')


# Get resids

# For every group, get 10k residuals, and take the average

fitted.vals = NULL
# For each group
for (j in 1:N) {
	index = ((j-1) * t + 1) : (t * j)
	
	betas = t(cbind(store.intercept[, j], store.price[, j], store.rainfall[, j]))
	X.j.B.j = apply(X[index, ] %*% betas, 1, mean)
	fitted.vals = c(fitted.vals, X.j.B.j)
}
resids = y - fitted.vals
