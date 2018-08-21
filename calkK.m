% stiffnesscalculator

function k = calcK (f_stopband,m)

k = (f_stopband*2*pi())^2*m
k1 = k/2
k2 = k/2
fe_SpringMass = sqrt((k1+k2)/m)/(2*pi())


end