using PyCall
@pyimport matplotlib.pyplot as plt

# Initial conditions
function from_prior(num::Int64)
	X = rand(num)
	Y = rand(num)
	return (X, Y)
end

# Upper corner count
function ucc(x::Float64, y::Float64, X::Vector{Float64}, Y::Vector{Float64})
	@assert length(X) == length(Y)
	return nothing
end

(X, Y) = from_prior(1000)

plt.plot(X, Y, "bo")
plt.axis("scaled")
plt.show()

