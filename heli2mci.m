function X_mci = heli2mci(X_heli)

% Input:
% X_heli: (3,1) vector in heliocentric coordinates
%
% Output:
% X_mci: (3,1) vector in MCI coordinates

alpha = 82.9076;
beta = 26.7180;
gamma = -177.6170;

C_alpha = [cosd(alpha) sind(alpha) 0;
           -sind(alpha) cosd(alpha) 0;
           0 0 1];

C_beta = [1 0 0;
          0 cosd(beta) sind(beta);
          0 -sind(beta) cosd(beta)];

C_gamma = [cosd(gamma) sind(gamma) 0;
           -sind(gamma) cosd(gamma) 0;
           0 0 1];


X_mci = C_gamma*C_beta*C_alpha*X_heli;

end