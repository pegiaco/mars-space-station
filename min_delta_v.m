function [min_delta_v_total, min_delta_v_date_launch, min_delta_v_date_arrival] = min_delta_v(N_launch, peak_indices, delta_v_total, date_launch_list, date_arrival_list)

% Input:
% N_launch: number of launch periods
% peak_indices: vector containing the indices of the peaks in Earth/Mars
% distance
% delta_v_total: matrix containing total Delta V
%
% Output:
% min_delta_v_total: vector containing the minimum total Delta V for each
% synodic period
% min_delta_v_date_launch: vector containing the launch date corresponding
% to the minimum total Delta V for each synodic period
% min_delta_v_date_arrival: vector containing the arrival date
% corresponding to the minimum total Delta V for each synodic period
% date_launch_list: datetime vector of launch dates
% date_arrival_list: datetime vector of arrival dates

min_delta_v_total = zeros(1, N_launch);
min_delta_v_date_launch = zeros(1, N_launch);
min_delta_v_date_arrival = zeros(1, N_launch);

for k = 1:N_launch
    % Indices corresponding to the current chunk
    chunck_indices = peak_indices(k):peak_indices(k+1)-1;

    % Extract the subset of delta_v_total for the current chunk
    delta_v_total_chunk = delta_v_total(chunck_indices, chunck_indices);

    % Find the minimum deltaV_total and its indices in the current chunk
    [min_delta_v, idx_min] = min(delta_v_total_chunk(:));
    [min_row, min_col] = ind2sub(size(delta_v_total_chunk), idx_min);

    % Save the associated launch date and arrival date
    min_launch_date = date_launch_list(chunck_indices(min_row));
    min_arrival_date = date_arrival_list(chunck_indices(min_col));

    % Save the results
    min_delta_v_total(k) = min_delta_v;
    min_delta_v_date_launch(k) = min_launch_date;
    min_delta_v_date_arrival(k) = min_arrival_date;
end

end