model = "model
{
  # Coefficients
  A ~ dunif(-1000, 1000)
  B ~ dunif(-1000, 1000)
  C ~ dunif(-1000, 1000)

  bias0 ~ dunif(-1000, 1000)
  slope ~ dunif(-10, 0)

  for(i in 1:N)
  {
    bias[i] <- bias0*walkers[i]^slope
    log_sigma[i] <- A*log(walkers[i]) + B*log(reps[i]) + C
    sigma[i] <- exp(log_sigma[i])
    # log(Z) as output variable
    y[i] ~ dnorm(-209.053 + bias[i], 1/sigma[i]^2)
    # H as output variable
    #y[i] ~ dnorm(48.4834, 1/sigma[i]^2)

    # Duplicate data for posterior predictive check
    y_copy[i] ~ dnorm(-209.053 + bias[i], 1/sigma[i]^2)
  }

  # Predict the results of an accurate run
  log_sigma_new <- A*log(100) + B*log(1000) + C
  sigma_new <- exp(log_sigma_new)
  y_new ~ dnorm(-209.053, 1/sigma_new^2)
}
"

data = as.matrix(read.table('results.txt'))

# Cull runs with < 3 walkers or reps
data = data[data[,1] >= 3 & data[,2] >= 3, ]

# log(Z) as output variable
data = list(walkers=data[,1], reps=data[,2], y=data[,4], N=length(data[,4]))
# H as output variable
# data = list(walkers=data[,1], reps=data[,2], y=data[,5], N=length(data[,4]))

# Variables to monitor
variable_names = c('A', 'B', 'C', 'bias0', 'slope', 'y_new', 'y_copy')

# How many burn-in steps?
burn_in = 1000

# How many proper steps?
steps = 100000

# Thinning?
thin = 10

# Random number seed
seed = 42


# NO NEED TO EDIT PAST HERE!!!
# Just run it all and use the results list.

library('rjags')

# Write model out to file
fileConn=file("model.temp")
writeLines(model, fileConn)
close(fileConn)

if(all(is.na(data)))
{
	m = jags.model(file="model.temp", inits=list(.RNG.seed=seed, .RNG.name="base::Mersenne-Twister"))
} else
{
	m = jags.model(file="model.temp", data=data, inits=list(.RNG.seed=seed, .RNG.name="base::Mersenne-Twister"))
}
update(m, burn_in)
draw = jags.samples(m, steps, thin=thin, variable.names = variable_names)
# Convert to a list
make_list <- function(draw)
{
	results = list()
	for(name in names(draw))
	{
		# Extract "chain 1"
		results[[name]] = as.array(draw[[name]][,,1])
		
		# Transpose 2D arrays
		if(length(dim(results[[name]])) == 2)
			results[[name]] = t(results[[name]])
	}
	return(results)
}
results = make_list(draw)

# Plot for posterior predictive check
plot(results$y_copy[100,], type='l', col='red', ylim=c(-220, -190))
lines(results$y_copy[200,], col='red')
lines(results$y_copy[300,], col='red')
lines(results$y_copy[400,], col='red')
lines(results$y_copy[500,], col='red')
lines(results$y_copy[600,], col='red')
lines(results$y_copy[700,], col='red')
lines(results$y_copy[800,], col='red')
lines(results$y_copy[900,], col='red')
lines(data$y, col='black', lwd=5)
