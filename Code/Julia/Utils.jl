# Heavy tailed proposal
function randh(size::Int64)
	return 10.0.^(1.5 - 6.0*rand(size)).*randn(size)
end

function logsumexp(x::Array{Float64, 1})
	top = maximum(x)
	y = exp(x - top)
	return log(sum(y)) + top
end

function logdiffexp(a::Float64, b::Float64)
	return b + log(exp(a - b) - 1.)
end

