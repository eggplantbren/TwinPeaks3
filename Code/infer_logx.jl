#using PyCall
#@pyimport matplotlib.pyplot as plt

# Load the files
scalars = transpose(readdlm("scalars.txt")) # Transpose for faster access
logw = vec(readdlm("logw.txt"))

# Truncate if necessary
N = minimum([size(scalars)[2], size(logw)[1]])
scalars = scalars[:, 1:N]
logw = logw[1:N]

# Which run each distribution belongs to
run_id = ones(Int64, N)
for i in 2:N
  if logw[i] == -1.
    run_id[i] = run_id[i-1] + 1
  else
    run_id[i] = run_id[i-1]
  end
end

# Generate logX values from the prior (i.e. not taking into account any
# inter-run comparisons)
logX = zeros(N)
for i in 1:N
  if logw[i] == -1.
    logX[i] = log(rand())
  else
    logX[i] = logX[i-1] + log(rand())
  end
end

# Gibbs sample one value
function update!(logx::Array{Float64, 1}, ii::Int64)
  # Generate proposal
  proposal = 0.

  # Range of values for proposal
  lower = -1E300
  upper = 0.

  # Lower limit -- next point in same run (if it exists)
  if (ii != N) && (run_id[ii+1] == run_id[ii])
    lower = logx[ii+1]
  end
  # Upper limit -- previous point in same run (if it exists)
  if (ii != 1) && (run_id[ii-1] == run_id[ii])
    upper = logx[ii-1]
  end

  # If lower limit exists, generate uniformly between limits
  if lower != -1E300
    proposal = lower + (upper - lower)*rand()
  else
    # Otherwise, use exponential distribution
    proposal = upper + log(rand())
  end

  # Measure inconsistency
  inconsistency_old = 0
  inconsistency_new = 0
  for j in 1:N
    if run_id[j] != run_id[ii]

      # Distribution ii is within distribution j
      if (logx[ii] > logx[j]) && all(scalars[:, ii] .>= scalars[:, j])
          inconsistency_old += 1
      end

      # Distribution ii is outside distribution j
      if (logx[ii] < logx[j]) && all(scalars[:, ii] .<= scalars[:, j])
	inconsistency_old += 1
      end

      # Distribution ii is within distribution j
      if (proposal > logx[j]) && all(scalars[:, ii] .>= scalars[:, j])
          inconsistency_new += 1
      end

      # Distribution ii is outside distribution j
      if (proposal < logx[j]) && all(scalars[:, ii] .<= scalars[:, j])
          inconsistency_new += 1
      end
    end
  end

  inconsistency = inconsistency_old
  # Accept?
  if inconsistency_new <= inconsistency_old
    logx[ii] = proposal
    inconsistency = inconsistency_new
  end
  return inconsistency
end

function update_all!(logx::Array{Float64, 1})
  total_inconsistency = 0
  for i in 1:N
    total_inconsistency += update!(logx, i)
  end
  return total_inconsistency
end

for i in 1:100
  inconsistency = update_all!(logX)
  println(i, " ", inconsistency)
end

writedlm("logX.txt", logX)

