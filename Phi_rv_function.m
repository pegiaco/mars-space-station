function Phi_rv = Phi_rv_function(n,t)

Phi_rv = zeros(2,2);

Phi_rv(1,1) = (1/n)*sin(n*t);
Phi_rv(1,2) = (2/n)*(1 - cos(n*t));
Phi_rv(2,1) = -Phi_rv(1,2);
Phi_rv(2,2) = (1/n)*(4*sin(n*t) - 3*n*t);

end