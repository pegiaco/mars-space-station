function [r_eci, v_eci] = ecef2eci(JD, r_ecef, v_ecef)

% Inputs:
% JD: Julian date
% r_ecef: position in ECEF coordinates [km]
% v_ecef: velocity in ECEF coordinates [km/sec]
%
% Outputs:
% r_eci: position in ECI coordinates [km]
% v_eci: velocity in ECI coordinates [km/sec]

w_Earth = 0.0000729211585530; % rad/sec

theta = GMST(JD);

eci_C_ecef = [cosd(theta), -sind(theta), 0;
         sind(theta), cosd(theta), 0;
         0, 0, 1];

r_eci = eci_C_ecef * r_ecef;

v_eci = v_ecef + cross([0;0;w_Earth], r_eci);

end