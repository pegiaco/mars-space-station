function Phi_vv = Phi_vv_function(n,t)

Phi_vv = zeros(2,2);

Phi_vv(1,1) = cos(n*t);
Phi_vv(1,2) = 2*sin(n*t);
Phi_vv(2,1) = -Phi_vv(1,2);
Phi_vv(2,2) = 4*cos(n*t) - 3;

end