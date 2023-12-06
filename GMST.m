function theta = GMST(JD)

% GMST: Greenwich Mean Sidereal Time

% Inputs:
% JD: Julian date
%
% Outputs:
% theta: GMST [deg]

global julian_offset;

w_Earth = 0.0000729211585530; % rad/sec

t_ref = [2014, 01, 01, 00, 00, 00];
JD_ref = datenum(t_ref)+julian_offset;
theta_ref = 6.7045617; % sidereal hours -- GMST at t_ref

theta = wrapTo360(theta_ref * 15 + (180/pi)*w_Earth*(JD-JD_ref)*86400);

end