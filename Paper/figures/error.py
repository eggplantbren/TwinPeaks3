from pylab import *

# Change font
rc("font", size=18, family="serif", serif="Computer Sans")
rc("text", usetex=True)

# Load results
results = loadtxt("../../Code/C++/results.dat")

# Plot first set of results
which = (results[:,1] == 500)
loglog(results[which, 0], results[which, 2], "bo-", label="500 MCMC steps")

## Plot second set of results
#which = (results[:,1] == 1000)
#loglog(results[which, 0], results[which, 2], "go-", label="1000 MCMC steps")

## Plot third set of results
#which = (results[:,1] == 2000)
#loglog(results[which, 0], results[which, 2], "yo-", label="2000 MCMC steps")

# Axis labels and legend
xlabel("Number of particles, $N$")
ylabel("RMS error in $\\ln Z(\\beta_1, \\beta_2)$")
#legend(loc="upper right")

# Save and display the figure
savefig("error.pdf", bbox_inches="tight")
show()

