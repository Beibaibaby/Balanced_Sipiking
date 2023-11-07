using Distributions
lambda_noise = 0.02
dt = 0.1

for i in 1:1000
  a=rand(Poisson(lambda_noise * dt))
  if a == 1
  println(a)
  end

end

println(Poisson(lambda_noise * dt))