include("Sampler.jl")

sampler = Sampler(100)
initialise!(sampler)

#while(true)
	do_iteration!(sampler)
#end

