using PyCall
@pyimport matplotlib.pyplot as plt
include("Utils.jl")

run(`./infer_logx`)

# Temperatures
T1, T2 = 0.1, 1.

# First calculate things about the scalars (e.g. the normalising constant)
scalars = readdlm("scalars.txt")
logw = vec(readdlm("logX.txt"))

smallest = minimum([size(scalars)[1], size(logw)[1]])
scalars = scalars[1:smallest, :]
logw = logw[1:smallest]

# Prior weights, normalised
logw = logw - logsumexp(logw)

# Posterior weights, unnormalised
logW = logw + scalars[:,1]./T1 + scalars[:,2]./T2

# Normaliser
logZ = logsumexp(logW)

# Posterior weights, normalised
logWW = logW - logZ
ess = exp(-sum(exp(logWW).*logWW))

# Information
H = sum(exp(logWW).*(logWW - logw))

println("log(Z) = ", logZ)
println("H = ", H, " nats.")

plt.plot(exp(logWW))
plt.ylabel("Weight wrt canonical distribution")
plt.title(string("ESS = ", ess))
plt.show()

