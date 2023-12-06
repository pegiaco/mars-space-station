function [] = porkchop_plot(n, title_str, filename, ...
                            date_launch_list, date_arrival_list, M, levels, ...
                            enable_chunks, peak_dates, margin, ...
                            enable_min, X_min, Y_min)

global output_folder;

figure(n);
% Plot the contour
hold on;
contour(date_launch_list, date_arrival_list, M', levels);
colorbar;
if enable_chunks == true
    % Plot the limits
    for i = 1:length(peak_dates)
        line([peak_dates(i), peak_dates(i)], ...
             get(gca, 'Ylim'), 'Color', 'red', 'LineStyle', '--');
        line(get(gca, 'Xlim'), [peak_dates(i)+margin, peak_dates(i)+margin], 'Color', 'blue', 'LineStyle', '--');
    end
end
if enable_min == true
    % Plot the location of the minimum values
    for i = 1:length(X_min)
        scatter(X_min(i), Y_min(i), 'filled', 'MarkerFaceColor', 'blue');
    end
end
hold off;
title(title_str + " - Earth to Mars mission planning");
xlabel("Earth departure (Julian date)");
ylabel("Mars arrival (Julian date)");
xlim([date_launch_list(1), date_launch_list(end)]);
ylim([date_arrival_list(1), date_arrival_list(end)]);
datetick('x', 2, 'keeplimits');
datetick('y', 2, 'keeplimits');
set(gca, 'XMinorTick', 'on');

% Save figure as PNG in a subfolder of the current directory
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
saveas(gcf, fullfile(output_folder, [filename, '_plot.png']));

end