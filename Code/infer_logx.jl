using PyCall
@pyimport matplotlib.pyplot as plt

# Generate initial guesses for logX to initialise the MCMC.
function starting_point(N::Int64, run_id::Array{Int64, 1})
  print("Generating initial guess...")
  logX = zeros(N)

  # Assign first run
  logX[1] = log(rand())
  k = 2
  while run_id[k] == 1
    logX[k] = logX[k-1] + log(rand())
    k += 1
  end

  for i in k:N
    # Limits
    upper = 0.
    lower = -1E300
    if run_id[i-1] == run_id[i]
      upper = logX[i-1]
    end

    for j in 1:(i-1)
      if run_id[j] != run_id[i]
        # If this other distribution is exterior, can reduce 'upper'
        if all(scalars[:, i] .>= scalars[:, j]) && logX[j] < upper
          upper = logX[j]
        end

        # If this other distribution is interior, can increase 'lower'
        if all(scalars[:, i] .<= scalars[:, j]) && logX[j] > lower
          lower = logX[j]
        end
      end
    end

    if lower == -1E300
      logX[i] = upper + log(rand())
    else
      logX[i] = lower + (upper - lower)*rand()
    end
  end
  println("done.")
  return logX
end

# Gibbs sample one value
function update!(logx::Array{Float64, 1}, ii::Int64, N::Int64,
                                scalars::Array{Float64, 2},
                                run_id::Array{Int64, 1})
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

function update_all!(logx::Array{Float64, 1}, N::Int64,
                                scalars::Array{Float64, 2},
                                run_id::Array{Int64, 1})
  total_inconsistency = 0
  for i in 1:N
    total_inconsistency += update!(logx, i, N, scalars, run_id)
    print('.')
  end
  return total_inconsistency
end

###################################################
# MAIN CODE
###################################################

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

logX = starting_point(N, run_id)

for i in 1:1000000
  inconsistency = update_all!(logX, N, scalars, run_id)
  println(i, " ", inconsistency)
  writedlm("logX.txt", logX)
end

