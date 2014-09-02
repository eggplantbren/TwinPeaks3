model = "model
{
  # Coefficients
  A ~ dunif(-1000, 1000)
  B ~ dunif(-1000, 1000)
  C ~ dunif(-1000, 1000)

  for(i in 1:N)
  {
    log_sigma[i] <- A*log(walkers[i]) + B*log(reps[i]) + C
    sigma[i] <- exp(log_sigma[i])
    # log(Z) as output variable
    y[i] ~ dnorm(-209.053, 1/sigma[i]^2)
    # H as output variable
    #y[i] ~ dnorm(48.4834, 1/sigma[i]^2)
  }
}
"

data = as.matrix(read.table('results.txt'))
# log(Z) as output variable
data = list(walkers=data[,1], reps=data[,2], y=data[,4], N=length(data[,4]))
# H as output variable
# data = list(walkers=data[,1], reps=data[,2], y=data[,5], N=length(data[,4]))

# Variables to monitor
variable_names = c('A', 'B', 'C')

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
