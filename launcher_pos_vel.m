function X = launcher_pos_vel(t, lambda, phi_gd, h, R, w, theta_g, t_0)

% position in ECEF cooridnates
r_ECEF = [(R+h)*cosd(phi_gd)*cosd(lambda);...
            (R+h)*cosd(phi_gd)*sind(lambda);...
            (R+h)*sind(phi_gd)];

% theta_g
theta_g = theta_g + (w*180/pi)*(t - t_0); % deg 

% position in ECI cooridnates
C = transpose([cosd(theta_g), sind(theta_g), 0;...
                -sind(theta_g), cosd(theta_g), 0;...
                0, 0, 1]);

r_ECI = C*r_ECEF;

% velocity in ECI coordinates
v_ECI = cross([0;0;w], r_ECI);


X = [r_ECI;v_ECI];

end