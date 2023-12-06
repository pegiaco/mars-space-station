function delta_V_f = delta_V_f_rendez_vous_fctn(n,t_f,X_0)

x_0 = X_0(1,1); %km
y_0 = X_0(2,1); %km
z_0 = X_0(3,1); %km

Phi_rv = Phi_rv_function(n,t_f);
Phi_rv_inv = inv(Phi_rv);
Phi_rr = Phi_rr_function(n,t_f);

Phi_vr = Phi_vr_function(n,t_f);
Phi_vv = Phi_vv_function(n,t_f);

temp = (Phi_vr - Phi_vv*Phi_rv_inv*Phi_rr)*[x_0; y_0];
x_dot_f_minus = temp(1,1);
y_dot_f_minus = temp(2,1);

z_dot_f_minus = -(n/sin(n*t_f))*z_0;

delta_V_f = [-x_dot_f_minus; -y_dot_f_minus; -z_dot_f_minus];

end