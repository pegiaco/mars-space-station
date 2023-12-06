function X_eci = heli2eci(X_heli)

% Input:
% X_heli: (3,1) vector in heliocentric coordinates
%
% Output:
% X_eci: (3,1) vector in ECI coordinates

epsilon = 23.45; % deg

eci_C_heli = [1, 0, 0; ...
              0, cosd(epsilon), -sind(epsilon); ...
              0, sind(epsilon), cosd(epsilon)];

X_eci = eci_C_heli * X_heli;

end