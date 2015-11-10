using PyCall
@pyimport matplotlib.pyplot as plt

@doc """
Initial conditions
""" ->
function from_prior(num::Int64)
	@assert num > 1
	x = rand(num)
	y = x.^2 + 0.1*randn(num)
	return hcat(x, y)
end

@doc """
Measure upper corner count of a single point
""" ->
function ucc(x::Float64, y::Float64, X::Vector{Float64}, Y::Vector{Float64})
	@assert length(X) == length(Y)
	count = 0
	for(i in 1:length(X))
		if((X[i] > x) && (Y[i] > y))
			count += 1
		end
	end
	return count
end

@doc """
Measure upper corner counts of a bunch of points
""" ->
function uccs(X::Vector{Float64}, Y::Vector{Float64})
	@assert length(X) == length(Y)
	counts = Array(Int64, (length(X), ))
	for(i in 1:length(X))
		counts[i] = ucc(X[i], Y[i], X, Y)
	end
	return counts
end

# Generate some points
params = from_prior(1000)

# Assign x and y
xy = -log(rand(size(params)))

# Plot
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(params[:,1], params[:,2], "b.")

plt.subplot(1, 2, 2)
plt.plot(xy[:,1], xy[:,2], "b.")
plt.show()


