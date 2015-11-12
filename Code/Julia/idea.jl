using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport colormaps
include("Utils.jl")

# Number of pixels
(ni, nj) = (1024, 1024)

# Prior distribution
log_prior = Array(Float64, (ni, nj))
for(j in 1:nj)
	x = (j - 0.5)/nj
	for(i in 1:ni)
		y = 1.0 - (i - 0.5)/ni
		log_prior[i, j] = -0.5*(x - 0.1)^2/0.1^2 - 0.5*(y - x^2)^2/0.1^2
	end
end
log_prior -= logsumexp(vec(log_prior))

# Upper corner mass
logG = Array(Float64, (ni, nj))
for(j in nj:-1:1)
	for(i in 1:ni)
		up = -Inf
		right = -Inf
		up_and_right = -Inf
		if(i != 1)
			up = logG[i-1, j]
		end
		if(j != nj)
			right = logG[i, j+1]
		end
		if((i != 1) && (j != nj))
			up_and_right = logG[i-1, j+1]
		end
		logG[i, j] = logsumexp([log_prior[i, j], up, right])
		logG[i, j] = logdiffexp(logG[i, j], up_and_right)
	end
end

plt.figure(figsize=(14, 6))
plt.subplot(1, 2, 1)
plt.imshow(exp(log_prior), cmap=colormaps.viridis, interpolation="nearest")
plt.title("\$\\pi(L_1, L_2)\$", fontsize=16)
plt.xlabel("\$L_1\$", fontsize=16)
plt.ylabel("\$L_2\$", fontsize=16)
plt.gca()[:set_xticklabels]([])
plt.gca()[:set_yticklabels]([])

plt.subplot(1, 2, 2)
plt.imshow(log_prior - logG, cmap=colormaps.viridis, interpolation="nearest")
plt.title("\$\\log\\left[\\pi(L_1, L_2)/G(L_1, L_2)\\right]\$", fontsize=16)
plt.xlabel("\$L_1\$", fontsize=16)
plt.ylabel("\$L_2\$", fontsize=16)
plt.gca()[:set_xticklabels]([])
plt.gca()[:set_yticklabels]([])
plt.show()

