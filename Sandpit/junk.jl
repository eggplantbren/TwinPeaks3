using PyCall
@pyimport matplotlib.pyplot as plt

@doc """
Initial conditions
""" ->
function from_prior(num::Int64)
	@assert num > 1
	X = rand(num)
	Y = rand(num)
	return (X, Y)
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
(X, Y) = from_prior(1000)

# Plot
plt.scatter(X, Y, marker="o", s=uccs(X, Y), alpha=0.2)
plt.axis("scaled")
plt.show()

