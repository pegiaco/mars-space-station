function Phi_vr = Phi_vr_function(n,t)

Phi_vr = zeros(2,2);

Phi_vr(1,1) = 3*n*sin(n*t);
Phi_vr(1,2) = 0;
Phi_vr(2,1) = 6*n*(cos(n*t) - 1);
Phi_vr(2,2) = 0;

end