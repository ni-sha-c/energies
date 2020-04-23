using ForwardDiff 
f(x::Vector) = sum(sin, x) + prod(tan, x)*sum(sqrt, x);
x = rand(5)
g = x -> ForwardDiff.gradient(f, x);
g(x)
