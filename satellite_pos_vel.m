function X = satellite_pos_vel(dt, M_0, n, e, mu, Omega, i, w)

% Given t, find M
M = M_0*pi/180 + n*dt;

% Find E by solving Kepler's equation
E = pi*1.5;
c = Inf; % convergence criterion
while c > 1e-3
    E_next = E + (M - E + e*sin(E))/(1 - e*cos(E));
    c = abs(E_next - E);
    E = E_next;
end

% Find nu
nu = 2*atan( sqrt((1+e)/(1-e)) * tan(E/2) );

% periphocal coordinates
a = (mu/(n^2))^(1/3);
p = a*(1 - e^2);

r_peri = [p*cos(nu)/(1+e*cos(nu)); p*sin(nu)/(1+e*cos(nu)); 0];
v_peri = [-sqrt(mu/p)*sin(nu); sqrt(mu/p)*(e+cos(nu)); 0];

% Convert to MCI coordinates
Omega = Omega*pi/180; % rad
i = i*pi/180; % rad
w = w*pi/180; % rad

C_1 = [cos(Omega), sin(Omega), 0;...
        -sin(Omega), cos(Omega), 0;...
        0, 0, 1];

C_2 = [1, 0, 0;...
        0, cos(i), sin(i);...
        0, -sin(i), cos(i)];

C_3 = [cos(w), sin(w), 0;...
        -sin(w), cos(w), 0;...
        0, 0, 1];

C = transpose(C_3*C_2*C_1);

r_MCI = C*r_peri;
v_MCI = C*v_peri;


X = [r_MCI;v_MCI];

end