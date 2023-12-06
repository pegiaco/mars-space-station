function delta_V_0 = delta_V_0_intercept_fctn(n,t_f,X_0)

x_0 = X_0(1,1); %km
y_0 = X_0(2,1); %km
z_0 = X_0(3,1); %km
x_dot_zero_minus = X_0(4,1); %km/s
y_dot_zero_minus = X_0(5,1); %km/s
z_dot_zero_minus = X_0(6,1); %km/s

Phi_rv = Phi_rv_function(n,t_f);
Phi_rv_inv = inv(Phi_rv);
Phi_rr = Phi_rr_function(n,t_f);

temp = -Phi_rv_inv*Phi_rr*[x_0; y_0];
x_dot_zero_plus = temp(1,1);
y_dot_zero_plus = temp(2,1);

z_dot_zero_plus = -(n*cos(n*t_f)/sin(n*t_f))*z_0;

delta_V_0 = [x_dot_zero_plus; y_dot_zero_plus; z_dot_zero_plus] - [x_dot_zero_minus; y_dot_zero_minus; z_dot_zero_minus];

end