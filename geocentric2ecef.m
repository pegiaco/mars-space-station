function [r_ecef, v_ecef] = geocentric2ecef(lambda, phi, h)

% Inputs:
% lambda: longitude [deg]
% phi: latitude [deg]
% h: elevation above Earth's sphere [km]
%
% Outputs:
% r_ecef: position in ECEF coordinates [km]
% v_ecef: velocity in ECEF coordinates [km/sec]

r_Earth = 6378.137; % km

r_ecef = [(r_Earth+h)*cosd(phi)*cosd(lambda);
          (r_Earth+h)*cosd(phi)*sind(lambda);
          (r_Earth+h)*sind(phi)];

v_ecef = [0; 0; 0];

end