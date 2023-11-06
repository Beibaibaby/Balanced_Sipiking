using Distributions
lambda_noise = 2
dt = 0.1

for i in 1:100
  println(rand(Poisson(lambda_noise * dt)))
end

println(Poisson(lambda_noise * dt))