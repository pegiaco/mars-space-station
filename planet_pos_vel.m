function X = planet_pos_vel(planet_id, julian_date, mu)

% Get the planetary elements at J2000
[elements, rates] = AA279j2000_planetary_elements(planet_id);

% Compute the planetary elements at the given Julian date

N = (julian_date - 2451545)/36525; % number of centuries since J2000

a = elements(1) + rates(1)*N; % AU
a = a*149597870.7; % km
e = elements(2) + rates(2)*N;
i = elements(3) + rates(3)*N; % deg
Omega = elements(4) + rates(4)*N; % deg
omega_bar = elements(5) + rates(5)*N; % deg
L = elements(6) + rates(6)*N; % deg

% Get M and w (omega)
M = L - omega_bar; % deg
w = omega_bar - Omega; % deg

% Convert angles in radians
i = i*pi/180; % rad
Omega = Omega*pi/180; % rad
M = M*pi/180; % rad
w = w*pi/180; % rad

% Find E by solving Kepler's equation
E = pi;
c = Inf; % convergence criterion
while c > 1e-3
    E_next = E + (M - E + e*sin(E))/(1 - e*cos(E));
    c = abs(E_next - E);
    E = E_next;
end

% Find nu
nu = 2*atan( sqrt((1+e)/(1-e)) * tan(E/2) );

% periphocal coordinates
p = a*(1 - e^2);

r_peri = [p*cos(nu)/(1+e*cos(nu)); p*sin(nu)/(1+e*cos(nu)); 0];
v_peri = [-sqrt(mu/p)*sin(nu); sqrt(mu/p)*(e+cos(nu)); 0];

% Convert to heliocentric coordinates

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

r_helio = C*r_peri;
v_helio = C*v_peri;


X = [r_helio;v_helio];

end