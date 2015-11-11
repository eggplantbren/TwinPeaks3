using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport colormaps

# Number of pixels
(ni, nj) = (101, 101)

# Prior distribution
prior = Array(Float64, (ni, nj))
for(j in 1:nj)
	for(i in 1:ni)
		prior[i, j] = exp(-0.5*(j - 0.5*nj)^2/10^2
							-0.5*(i - 0.5*ni)^2/10^2)
	end
end
prior = prior/sum(prior)

# Upper corner mass
G = Array(Float64, (ni, nj))
for(j in nj:-1:1)
	for(i in 1:ni)
		up = 0.0
		right = 0.0
		up_and_right = 0.0
		if(i != 1)
			up += G[i-1, j]
		end
		if(j != nj)
			right += G[i, j+1]
		end
		if((i != 1) && (j != nj))
			up_and_right += G[i-1, j+1]
		end

		G[i, j] = prior[i, j] + up + right - up_and_right
	end
end

plt.figure(figsize=(14, 6))
plt.subplot(1, 2, 1)
plt.imshow(prior, cmap=colormaps.viridis, interpolation="nearest")
plt.title("\$\\pi(L_1, L_2)\$")
plt.xlabel("\$L_1\$")
plt.ylabel("\$L_2\$")
plt.gca()[:set_xticklabels]([])
plt.gca()[:set_yticklabels]([])

plt.subplot(1, 2, 2)
plt.imshow(prior./G, cmap=colormaps.viridis, interpolation="nearest")
plt.title("\$\\pi(L_1, L_2)/G(L_1, L_2)\$")
plt.xlabel("\$L_1\$")
plt.ylabel("\$L_2\$")
plt.gca()[:set_xticklabels]([])
plt.gca()[:set_yticklabels]([])


plt.show()

