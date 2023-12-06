function Phi_rr = Phi_rr_function(n,t)

Phi_rr = zeros(2,2);

Phi_rr(1,1) = 4 - 3*cos(n*t);
Phi_rr(1,2) = 0;
Phi_rr(2,1) = 6*(sin(n*t) - n*t);
Phi_rr(2,2) = 1;

end