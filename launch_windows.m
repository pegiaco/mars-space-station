clear all; close all; clc;

%% Constants

global mu_Sun;
mu_Sun = 1.3271244004193938e11; % km^3/sec^2

global julian_offset;
julian_offset = 1721058.5;

mu_Earth = 398600; % km^3/sec^2
global mu_Mars;
mu_Mars = 42828; % km^3/sec^2
r_Earth_Mars = 225e6; % km -- Mean Earth-Mars distance

% Synodic period
T_Earth = 365.25; % days/rev
T_Mars = 687; % days/rev
T_syn = 1/(1/T_Earth - 1/T_Mars); % days/rev

margin = T_syn/3;

% Launch site (KOUROU) - Geocentric coordinates
lambda = -52.650299; % deg
phi = 5.159700; % deg
h = -14.7e-3; % km

R_Earth = 6378; % km
R_Mars = 3389.5; % km

h_orbit_Earth = 400; % km
h_orbit_Mars = 1000; % km

% Subfolder for plots
global output_folder;
output_folder = 'plots';


%% Julian date ranges

N_launch = 3;

step = 5; % days

% Start and end dates for Earth departure
start_launch = [2024, 00, 00, 00, 00, 00];
Earth_start_date = datenum(start_launch) + julian_offset; % Julian date
Earth_end_date = datenum(start_launch) + (N_launch+1)*T_syn + T_syn/5 + julian_offset; % Julian date

% Date ranges
JD_Earth_departure = Earth_start_date:step:Earth_end_date;
JD_Mars_arrival = JD_Earth_departure + margin;


%% Compute Earth-Mars distance

Earth_Mars_distance = zeros(length(JD_Earth_departure), 1);

for i = 1:length(JD_Earth_departure)
    % Get planet positions
    Earth_pos_vel_temp = planet_pos_vel(3, JD_Earth_departure(i), mu_Sun);
    Mars_pos_vel_temp = planet_pos_vel(4, JD_Earth_departure(i), mu_Sun);

    % Earth-Mars distance
    Earth_Mars_distance(i) = norm(Mars_pos_vel_temp(1:3) - Earth_pos_vel_temp(1:3));
end

save('results/Earth_Mars_distance.mat', 'Earth_Mars_distance');

%% Earth-Mars distance plot

load('results/Earth_Mars_distance.mat');

date_launch_list = datenum(datetime(JD_Earth_departure, 'ConvertFrom', 'juliandate'));
date_arrival_list = datenum(datetime(JD_Mars_arrival, 'ConvertFrom', 'juliandate'));


% Find peaks in the Earth-Mars distance data
[peaks, peak_indices] = findpeaks(Earth_Mars_distance, 'MinPeakDistance', 10);

% Convert peak indices to corresponding dates
peak_dates = date_launch_list(peak_indices);

figure(1);
plot(date_launch_list, Earth_Mars_distance);
% Plot the peaks
hold on;
plot(peak_dates, peaks, 'ro', 'MarkerFaceColor', 'r');
for i = 1:length(peak_indices)
    line([peak_dates(i), peak_dates(i)], ...
         get(gca, 'Ylim'), 'Color', 'red', 'LineStyle', '--');
end
hold off;
title("Earth-Mars distance - Earth to Mars mission planning");
xlabel("Earth departure (Julian date)");
ylabel("Distance (km)");
xlim([date_launch_list(1), date_launch_list(end)]);
datetick('x', 2, 'keeplimits');
set(gca, 'XMinorTick', 'on');

% Save figure as PNG in a subfolder of the current directory
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
saveas(gcf, fullfile(output_folder, ['Earth_Mars_distance', '_plot.png']));


%% Compute data

% Initialize data structures
transfer_time = nan(length(JD_Earth_departure), length(JD_Mars_arrival));
C3 = nan(length(JD_Earth_departure), length(JD_Mars_arrival));
RLA = nan(length(JD_Earth_departure), length(JD_Mars_arrival));
DLA = nan(length(JD_Earth_departure), length(JD_Mars_arrival));
v_plus_P_Earth = nan(length(JD_Earth_departure),length(JD_Mars_arrival));
i_launch = nan(length(JD_Earth_departure), length(JD_Mars_arrival));

v_Mars_minus_inf_mci = nan(3, length(JD_Earth_departure), length(JD_Mars_arrival));
v_minus_P_Mars = nan(length(JD_Earth_departure),length(JD_Mars_arrival));
i_insertion_MCI = nan(length(JD_Earth_departure),length(JD_Mars_arrival));

delta_v_launch = nan(length(JD_Earth_departure),length(JD_Mars_arrival));
delta_v_circular = nan(length(JD_Earth_departure),length(JD_Mars_arrival));
delta_v_equatorial = nan(length(JD_Earth_departure),length(JD_Mars_arrival));
delta_v_total = nan(length(JD_Earth_departure),length(JD_Mars_arrival));

for k = 1:N_launch
    % Indices corresponding to the current chunk
    chunck_indices = peak_indices(k):peak_indices(k+1)-1;

    for i = chunck_indices
        for j = chunck_indices
            if JD_Mars_arrival(j) > JD_Earth_departure(i)
        
                % Check if within desired range
                if date_arrival_list(j) >= peak_dates(k) + margin && ...
                   date_arrival_list(j) < peak_dates(k+1) + margin
        
                    % Transfer time
                    transfer_time(i,j) = JD_Mars_arrival(j) - JD_Earth_departure(i); % days
                    dtsec = transfer_time(i,j)*86400; % sec
                    
                    % Get planet positions
                    Earth_pos_vel_temp = planet_pos_vel(3, JD_Earth_departure(i), mu_Sun);
                    Mars_pos_vel_temp = planet_pos_vel(4, JD_Mars_arrival(j), mu_Sun);

                    c = norm(Mars_pos_vel_temp(1:3) - Earth_pos_vel_temp(1:3)); % km
                    s = (norm(Earth_pos_vel_temp(1:3)) + norm(Mars_pos_vel_temp(1:3)) + c)/2; % km
                    t_parabola = (1/3)*sqrt(2/mu_Sun)*(s^(3/2) - (s-c)^(3/2)); % sec

                    % Check if parabola or hyperbola
                    if t_parabola < transfer_time(i,j)*86400
                        % Solve Lambert problem
                        [v_Earth_plus_temp, v_Mars_minus_temp, error] = AA279lambert_curtis(mu_Sun, ...
                            Earth_pos_vel_temp(1:3), Mars_pos_vel_temp(1:3), 'pro', 0, dtsec);
            
                        % Earth departure velocity in ECI coordinates
                        v_plus_inf_heli_temp = v_Earth_plus_temp - Earth_pos_vel_temp(4:6);
                        v_plus_inf_eci_temp = heli2eci(v_plus_inf_heli_temp);
            
                        % Hyperbolic launch parameters
                        C3(i,j) = norm(v_plus_inf_eci_temp)^2; % km^2/sec^2
                        RLA(i,j) = asind(v_plus_inf_eci_temp(3)/norm(v_plus_inf_eci_temp)); % deg
                        DLA(i,j) = atan2d(v_plus_inf_eci_temp(2), v_plus_inf_eci_temp(1)); % deg
                        v_plus_P_Earth(i,j) = sqrt(norm(v_plus_inf_heli_temp)^2+2*mu_Earth/(R_Earth+h_orbit_Earth));
                        theta_GMST = GMST(JD_Earth_departure(i));
                        alpha_L = lambda + theta_GMST;
                        % Select positive solution for beta (prograde/East launch)
                        beta = mod(atan2d(sind(RLA(i,j) - alpha_L), ...
                                          cosd(phi)*tand(DLA(i,j)) - sind(phi)*cosd(RLA(i,j) - alpha_L)), 360);
                        i_launch(i,j) = acosd(cosd(phi)*sind(beta));

                        % Mars arrival velocity
                        v_Mars_minus_inf_heli_temp = v_Mars_minus_temp - Mars_pos_vel_temp(4:6);
                        v_Mars_minus_inf_mci_temp = heli2mci(v_Mars_minus_inf_heli_temp);
                        v_Mars_minus_inf_mci(:,i,j) = v_Mars_minus_inf_mci_temp;

                        % Hyperbolic insertion parameters 
                        v_minus_P_Mars(i,j) = sqrt(norm(v_Mars_minus_inf_heli_temp)^2+2*mu_Mars/(R_Mars+h_orbit_Mars));
                        i_insertion_MCI(i,j) = asind(v_Mars_minus_inf_mci_temp(3)/norm(v_Mars_minus_inf_mci_temp));
                        
                        % Delta V
                        delta_v_launch(i,j) = abs(v_plus_P_Earth(i,j)-sqrt(mu_Earth/(R_Earth+h_orbit_Earth)));
                        delta_v_circular(i,j) = abs(sqrt(mu_Mars/(R_Mars+h_orbit_Mars))-v_minus_P_Mars(i,j));
                        delta_v_equatorial(i,j) = sqrt(mu_Mars/(R_Mars+h_orbit_Mars))*abs(pi/180*i_insertion_MCI(i,j));
                        delta_v_total(i,j) = delta_v_launch(i,j)+delta_v_circular(i,j)+delta_v_equatorial(i,j);

                    end
                end
            end
        end
    end
end

save('results/transfer_time.mat','transfer_time');
save('results/C3.mat','C3');
save('results/RLA.mat','RLA');
save('results/DLA.mat','DLA');
save('results/v_plus_P_Earth.mat','v_plus_P_Earth');
save('results/i_launch.mat','i_launch');
save('results/v_minus_P_Mars.mat','v_minus_P_Mars');
save('results/i_insertion_MCI.mat','i_insertion_MCI');
save('results/delta_v_launch.mat','delta_v_launch');
save('results/delta_v_circular.mat','delta_v_circular');
save('results/delta_v_equatorial.mat','delta_v_equatorial');
save('results/delta_v_total.mat','delta_v_total');


%% Entire Delta V total porkchop plot

load('results/delta_v_total.mat');

porkchop_plot(4+3*N_launch+1, "\DeltaV_{total}", 'Delta_V_total_all', ...
              date_launch_list, date_arrival_list, delta_v_total, 5:1:20, ...
              true, peak_dates, margin, ...
              false, [], []);


%% Individual Delta V total porkchop plots

load('results/delta_v_total.mat');

for k = 1:N_launch
    % Indices corresponding to the current chunk
    chunck_indices = peak_indices(k):peak_indices(k+1)-1;

    porkchop_plot(5+3*N_launch+k, "\DeltaV_{total} (window #"+k+")", ['Delta_V_total_window_', num2str(k)], ...
                  date_launch_list(chunck_indices), date_arrival_list(chunck_indices), delta_v_total(chunck_indices, chunck_indices), 5:1:20, ...
                  false, peak_dates, margin, ...
                  false, [], []);
end


%% Entire Delta V total porkchop plot with minimum values

load('results/delta_v_total.mat');

[min_delta_v_total, min_delta_v_date_launch, min_delta_v_date_arrival] = ...
    min_delta_v(N_launch, peak_indices, delta_v_total, date_launch_list, date_arrival_list);

porkchop_plot(5+4*N_launch+1, "\DeltaV_{total}", 'Delta_V_total_with_min_all', ...
              date_launch_list, date_arrival_list, delta_v_total, 5:1:20, ...
              true, peak_dates, margin, ...
              true, min_delta_v_date_launch, min_delta_v_date_arrival);


%% Individual Delta V total porkchop plots with minimum values

load('results/delta_v_total.mat');

[min_delta_v_total, min_delta_v_date_launch, min_delta_v_date_arrival] = ...
    min_delta_v(N_launch, peak_indices, delta_v_total, date_launch_list, date_arrival_list);

for k = 1:N_launch
    % Indices corresponding to the current chunk
    chunck_indices = peak_indices(k):peak_indices(k+1)-1;

    porkchop_plot(6+4*N_launch+k, "\DeltaV_{total} (window #"+k+")", ['Delta_V_total_with_min_window_', num2str(k)], ...
                  date_launch_list(chunck_indices), date_arrival_list(chunck_indices), delta_v_total(chunck_indices, chunck_indices), 5:1:20, ...
                  false, peak_dates, margin, ...
                  true, [min_delta_v_date_launch(k)], [min_delta_v_date_arrival(k)]);
end


%% Modules docking

% Load the Mars texture image
marsTexture = imread('mars.jpg');

% Mean motion
n_module = sqrt(mu_Mars/((R_Mars+h_orbit_Mars)^3)); % rad/sec

% Orbital period
T_module = 2*pi/n_module; % sec

% Initial orbit station
i_station = 0; % deg
Omega_station = 0; % deg
e_station = 0.2;
w_station = 40; % deg
M0_station = 0; % deg

a_station = 1.75*R_Mars; % km
T_station = 2*pi*sqrt((a_station^3)/mu_Mars); % sec
n_station = sqrt(mu_Mars/(a_station^3)); % sec


window = 17/24; % days
date_step = 1/(24*240); % days

t0_docking_sim = zeros(1, N_launch);

transfer_delta_v = zeros(1, N_launch);
transfer_start_dates = zeros(1, N_launch);
transfer_end_dates = zeros(1, N_launch);


for k=1:N_launch
    % Assumption: 1 entire orbit is traveled after the last Delta V for 0 deg inclination
    t0_docking_sim(k) = min_delta_v_date_arrival(k) + T_module/86400; % days

    start_dates = t0_docking_sim(k):date_step:(t0_docking_sim(k) + window); % days
    end_dates = t0_docking_sim(k):date_step:(t0_docking_sim(k) + window); % days

    v_module_start = nan(length(start_dates), length(end_dates), 3);
    delta_v_start = nan(length(start_dates), length(end_dates));
    delta_v_rdv = nan(length(start_dates), length(end_dates));
    delta_v_docking = nan(length(start_dates), length(end_dates));

    for i = 1:length(start_dates)
        % Position of the module in MCI coordinates
        t_module = (start_dates(i)-t0_docking_sim(k))*86400; % sec
        theta_module = n_module*t_module; % rad
        p_module = [(R_Mars+h_orbit_Mars)*cos(theta_module), (R_Mars+h_orbit_Mars)*sin(theta_module), 0]; % km
        v_module = [sqrt(mu_Mars/(R_Mars + h_orbit_Mars))*(-sin(theta_module)), sqrt(mu_Mars/(R_Mars + h_orbit_Mars))*cos(theta_module), 0]; % km/sec

        for j = 1:length(end_dates)
            % Assumption: a fictive station is at theta_station = 0 deg at t0_docking_sim(1)
            t_station = end_dates(j)*86400 - t0_docking_sim(1)*86400; % sec
            X = satellite_pos_vel(t_station, M0_station, n_station, e_station, mu_Mars, Omega_station, i_station, w_station);
            p_station = X(1:3);
            v_station = X(4:6);

            c = norm(p_station - p_module);
            s = (norm(p_module) + norm(p_station) + c)/2;
            
            dtsec = (end_dates(j) - start_dates(i))*86400; % sec

            % % Choose between short-way and long-way solutions
            % h_n = cross(p_module, v_module);
            % Tran_n = cross(p_module, p_station);
            % 
            % if dot(h_n, Tran_n) > 0
                % Short-way solution is better
                delta_t_para_short = (1/3)*sqrt(2/mu_Mars)*(s^(3/2) - (s-c)^(3/2)); % short-way

                % Check parabola or hyperbola (short-way)
                if dtsec > delta_t_para_short
                    % Solve Lambert problem (short-way solution)
                    [v1_short, v2_short, error] = AA279lambert_vallado_u(mu_Mars, ...
                                          p_module, p_station, 's', 0, dtsec);

                    v1_sim = v1_short;
                    v2_sim = v2_short;

                    % Check if the module is hitting Mars
                    if (dot(p_module,v1_sim) < 0 && dot(p_station,v2_sim) > 0)  % short-way solutions only
                        E_sim = (norm(v1_sim)^2)/2 - mu_Mars/norm(p_module);
                        a_sim = -mu_Mars/(2*E_sim);
                        h_sim = norm(cross(p_module,v1_sim));
                        p_sim = h_sim^2/mu_Mars;
                        e_sim = sqrt((a_sim-p_sim)/a_sim);
                        r_p = a_sim*(1-e_sim);
                        if r_p <= R_Mars % + h_orbit_limit
                            % non-physical values already set to NaN
                        else
                            % no collision
                            v_module_start(i,j,:) = v1_sim;
                            delta_v_start(i,j) = norm(v1_sim - v_module); % km/sec
                            delta_v_rdv(i,j) = norm(v_station - v2_sim); % km/sec
                            delta_v_docking(i,j) = delta_v_start(i,j) + delta_v_rdv(i,j); % sec
                        end
                    else
                        % no collision
                        v_module_start(i,j,:) = v1_sim;
                        delta_v_start(i,j) = norm(v1_sim - v_module); % km/sec
                        delta_v_rdv(i,j) = norm(v_station - v2_sim); % km/sec
                        delta_v_docking(i,j) = delta_v_start(i,j) + delta_v_rdv(i,j); % sec
                    end
                else
                   % non-physical values already set to NaN
                end

        end
    end

    % Find the minimum total Delta V for docking
    [min_delta_v, idx_min_delta_v] = min(delta_v_docking, [], 'all');
    [min_idx_i, min_idx_j] = ind2sub(size(delta_v_docking), idx_min_delta_v);

    % Extract the corresponding information
    min_delta_v_start = delta_v_start(min_idx_i, min_idx_j, :);
    min_docking_date_start = start_dates(min_idx_i); % days
    min_docking_date_rdv = end_dates(min_idx_j); % days

    transfer_delta_v(k) = min_delta_v_start;
    transfer_start_dates(k) = min_docking_date_start;
    transfer_end_dates(k) = min_docking_date_rdv;

    % Porkchop plot
    figure(200+k);
    set(gcf, 'Position', [400, 150, 900, 800]);
    contour(start_dates, end_dates, delta_v_docking', 0:0.1:20);
    colorbar;
    hold on;
    scatter(start_dates(min_idx_i), end_dates(min_idx_j), 'filled', 'red');
    hold off;
    title("\DeltaV_{total} - Docking (Module #"+k+")");
    xlabel('Start date (Julian date)');
    ylabel('End date (Julian date)');
    % Convert Julian dates to date strings
    dateStrings = datestr(start_dates);
    % Set x-axis ticks and labels
    xTicks = linspace(min(start_dates), max(start_dates), 6); % Adjust the number of ticks as needed
    xticks(xTicks);
    xticklabels(datestr(xTicks));
    % Set y-axis ticks and labels
    yTicks = linspace(min(end_dates), max(end_dates), 6); % Adjust the number of ticks as needed
    yticks(yTicks);
    yticklabels(datestr(yTicks));
    set(gca,'XMinorTick','on','YMinorTick','on')
    axis equal;
    saveas(gcf, fullfile(output_folder, ['Delta_v_total_docking_module_', num2str(k), '_plot.png']));


    
    % Position of the module at min_docking_date_start
    t_module = (start_dates(min_idx_i) - t0_docking_sim(k))*86400; % sec
    theta_module = n_module*t_module; % rad
    p_module = [(R_Mars+h_orbit_Mars)*cos(theta_module), (R_Mars+h_orbit_Mars)*sin(theta_module), 0]; % km

    % Trajectory of the module
    X_0 = [p_module(1); p_module(2); p_module(3); ...
           v_module_start(min_idx_i,min_idx_j,1);...
           v_module_start(min_idx_i,min_idx_j,2);...
           v_module_start(min_idx_i,min_idx_j,3)];
    dtsec = (end_dates(min_idx_j) - start_dates(min_idx_i))*86400; % sec
    transfer_step = 0.1; % sec
    t_span = 0:transfer_step:dtsec; % sec
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [t_out, X_out] = ode113(@state_derivative_2_body, t_span', X_0, options);

    % Needed to plot the ellipse
    T_ellipse = 0:1:T_station;
    x_station_ellipse = zeros(1, length(T_ellipse));
    y_station_ellipse = zeros(1, length(T_ellipse));
    z_station_ellipse = zeros(1, length(T_ellipse));
    for q = 1:length(T_ellipse)
        X = satellite_pos_vel(T_ellipse(q), M0_station, n_station, e_station, mu_Mars, Omega_station, i_station, w_station);
        x_station_ellipse(q) = X(1);
        y_station_ellipse(q) = X(2);
        z_station_ellipse(q) = X(3);
    end

    % Create a figure for the animation
    figure(300 + k);
    set(gcf, 'Position', [100, 150, 800, 800]);

    % Determine the time span for animation
    animationStartTime = t0_docking_sim(k);
    animationEndTime = end_dates(min_idx_j) + 2*T_station/86400; % You can adjust this based on your preference

    sec_per_frame = 30;
    numFrames = (animationEndTime-animationStartTime)*86400/sec_per_frame;

    % Loop over time steps and update the animation
    for frame = 0:1:(numFrames-1)
        % Plot the reference frame using quiver
        vector_length = 2*R_Mars;
        quiver3(0, 0, 0, vector_length, 0, 0, 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
        hold on;
        quiver3(0, 0, 0, 0, vector_length, 0, 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
        quiver3(0, 0, 0, 0, 0, vector_length, 'k', 'LineWidth', 1.5, 'HandleVisibility','off');

        % Plot Mars with texture
        [x, y, z] = sphere(50);
        surf(x * R_Mars, y * R_Mars, z * R_Mars, 'CData', marsTexture, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'HandleVisibility','off');

        % Plot the ellipse for the orbit of the station
        plot3(x_station_ellipse, y_station_ellipse, z_station_ellipse, 'Color', 'green', 'LineWidth', 0.5);

        % Plot the ellipse axis
        midpoint_idx = ceil(length(T_ellipse)/2);
        plot3(x_station_ellipse([1,midpoint_idx]), y_station_ellipse([1,midpoint_idx]), z_station_ellipse([1,midpoint_idx]), ...
            'Color', 'black', 'LineWidth', 0.2, 'LineStyle', '--', 'HandleVisibility','off');

        % Plot the circle for the orbit of the module
        theta_module_circle = linspace(0, 2*pi, 100);
        x_module_circle = (R_Mars + h_orbit_Mars) * cos(theta_module_circle);
        y_module_circle = (R_Mars + h_orbit_Mars) * sin(theta_module_circle);
        z_module_circle = zeros(size(theta_module_circle));
        plot3(x_module_circle, y_module_circle, z_module_circle, 'Color', 'red', 'LineWidth', 0.5);

        % Plot the trajectory of the module
        plot3(X_out(:,1), X_out(:,2), X_out(:,3), 'Color', 'blue', 'LineWidth', 1.1);

        % Calculate the animation time in days
        t_anim = animationStartTime*86400 + frame*sec_per_frame;

        % Plot the position of the station at the current time step
        t_station_anim = t_anim - t0_docking_sim(1) * 86400;
        X = satellite_pos_vel(t_station_anim, M0_station, n_station, e_station, mu_Mars, Omega_station, i_station, w_station);
        p_station_anim = X(1:3);
        if t_anim/86400 < end_dates(min_idx_j)
            scatter3(p_station_anim(1), p_station_anim(2), p_station_anim(3), 40, 'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
        else
            scatter3(p_station_anim(1), p_station_anim(2), p_station_anim(3), 40, 'filled', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k');
        end

        % Plot the position of the module at the current time step
        theta_module_anim = n_module * (frame*sec_per_frame);
        if t_anim < start_dates(min_idx_i)*86400
            p_module_anim = [(R_Mars + h_orbit_Mars) * cos(theta_module_anim), (R_Mars + h_orbit_Mars) * sin(theta_module_anim), 0];
            scatter3(p_module_anim(1), p_module_anim(2), p_module_anim(3), 40, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
        elseif t_anim > start_dates(min_idx_i)*86400 && t_anim < end_dates(min_idx_j)*86400
            p_module_anim = X_out(ceil((t_anim-start_dates(min_idx_i)*86400)/transfer_step), 1:3);
            scatter3(p_module_anim(1), p_module_anim(2), p_module_anim(3), 40, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
        else
            % Module docked to station, handled above
        end

        hold off;

        % Set axis properties
        axis equal;
        title(['Docking (Module #' num2str(k) ') - Date: ' datestr(t_anim/86400)]);
        xlabel('X (km)');
        ylabel('Y (km)');
        zlabel('Z (km)');
        xlim([-2.1, 2.1]*R_Mars);
        ylim([-2.1, 2.1]*R_Mars);
        % Set view from the top
        view(0, 90);
        % Add legend
        legend('Station orbit', 'Module orbit', 'Transfer orbit', 'Station', 'Module', 'Location', 'northeast');

        if frame == 0
            saveas(gcf, fullfile(output_folder, ['Transfer_module_', num2str(k), '_plot.png']));
        end

        % Capture the current frame for the video
        frames(frame+1) = getframe(gcf);

        % Remove the previous station position to update it in the next iteration
        clf;
    end

    % Create a video from the captured frames
    video = VideoWriter(['videos/docking_animation_', num2str(k), '.mp4'], 'MPEG-4');
    video.FrameRate = 60;

    open(video);
    writeVideo(video, frames);
    close(video);

end
