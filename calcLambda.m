function lambda = calcLambda(f,E,rho,h)
lambda = nthroot(E*h^2/(rho*12),4)*sqrt(2*pi/f)
end