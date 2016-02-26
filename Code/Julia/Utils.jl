@doc """
Heavy tailed distribution used for proposals
""" ->
function randh()
	return 10.0^(1.5 - 6.0*rand())*randn()
end

@doc """
log(sum(exp(x)))
""" ->
function logsumexp(x::Array{Float64, 1})
	top = maximum(x)
	y = exp(x - top)
	return log(sum(y)) + top
end

@doc """
log(exp(a) - exp(b))
""" ->
function logdiffexp(a::Float64, b::Float64)
	@assert a > b
	if(b == -Inf)
		return a
	end
	return b + log(exp(a - b) - 1.)
end

@doc """
Compare x and y, both vectors of floats (scalars).
Returns 1 if x is higher, -1 if x is lower, and 0 if they're not
comparable.
""" ->
function compare(x::Vector{Float64}, y::Vector{Float64})
    @assert length(x) == length(y)
    if all(x .> y)
        return 1
    elseif all(x .< y)
        return -1
    end
    return 0
end

