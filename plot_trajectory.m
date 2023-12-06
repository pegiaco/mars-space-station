function [] = plot_trajectory(n, t_f, X_0, delta_V_0, markerFrequency, figureNumber)

x_0 = X_0(1,1);
y_0 = X_0(2,1);

T = 0:1:t_f;
X = zeros(size(T));
Y = zeros(size(T));

for i = 1:length(T)
    temp = Phi_rr_function(n,T(i))*[x_0;y_0] + Phi_rv_function(n,T(i))*delta_V_0(1:2,1);
    X(i) = temp(1,1);
    Y(i) = temp(2,1);
end

figure(figureNumber);
plot(Y,X);
hold on;
% Chaser initial position
scatter(y_0, x_0, 'r', 'filled');
% Target position
scatter(0, 0, 'g', 'filled');
% Markers
for i = 1:length(T)
    if rem(T(i),markerFrequency*60) == 0
        scatter(Y(i),X(i),'blue');
        hold on;
    end
end
hold off;
title("Trajectory for rendez-vous after t_f = " + t_f/60 + " min | Markers every " + markerFrequency + " min");
xlabel("y");
ylabel("x");
axis equal;
grid on;
set(gca,'XDir','reverse');

end

