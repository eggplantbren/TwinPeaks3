using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport colormaps

# Number of pixels
(ni, nj) = (21, 21)

# Pixel widths
(dx, dy) = (1.0/nj, 1.0/ni)

# Prior distribution
prior = Array(Float64, (ni, nj))
for(j in 1:nj)
	for(i in 1:ni)
		prior[i, j] = exp(-0.5*(i - 0.5*ni)^2/2^2 - 0.5*(j - (ni - i))^2/2^2)
	end
end
prior = prior/sum(prior)

# Upper corner mass
G = Array(Float64, (ni, nj))
for(j in 1:nj)
	for(i in 1:ni)
		G[i, j] = 0.0
		for(jj in j:nj)
			for(ii in 1:i)
				G[i, j] += prior[ii, jj]
			end
		end
	end
end

plt.imshow((prior./G), cmap=colormaps.viridis, interpolation="nearest")
plt.show()

