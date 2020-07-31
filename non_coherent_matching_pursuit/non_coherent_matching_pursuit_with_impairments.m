clear;
run('data_for_compressive_sensing')

%% Basic parameters
alpha = linspace(-pi/2, pi/2, 1024);
% Vector gain_error_options and phase_error_options should be of the same length (just for simulation simplicity)
gain_error_options = 0:0.5:5; % In [dB]
phase_error_options = 0:5:50; % In [degree]
M_options = 5:5:15;
SNR = 20; % In [dB]


%% AoA estimation
filename = [save_dir,'measurement_raw_',num2str(SNR),'dB.csv'];
y = csvread(filename); % Read noisy data at specified SNR

AoA_RMSE = zeros(length(M_options), length(gain_error_options), 2); % Third dimension denotes the number of considered impairments (gain and phase)
detection_probability = zeros(length(M_options), length(gain_error_options), 2);

for error_index = 1:length(gain_error_options)
    gain_error_std = gain_error_options(error_index);
    phase_error_std = phase_error_options(error_index)/180*pi;
    
    %% Generate errors
    % Gain error
    gain_errors_dB = randn(Nr, max(M_options))*gain_error_std; % Gaussian distribution
    gain_errors_lin = 10.^(gain_errors_dB/10);
    % Phase error
    phase_errors = randn(Nr, max(M_options))*phase_error_std; % Gaussian distribution

    %% Consider different number of measurements
    for m_index = 1:length(M_options)
        M = M_options(m_index);

        %% Add hardware impairments separately: gain and phase error
        W_g = gain_errors_lin(:,1:M).*W(:,1:M);
        W_p = exp(1j*phase_errors(:,1:M)).*W(:,1:M);

        %% Generate non-coherent codebooks with impairments
        W_nc_g = zeros(M, length(alpha));
        W_nc_p = zeros(M, length(alpha));
        for i = 1:length(alpha)
            arx = exp(-1j*(0:Nr-1)'*pi*sin(alpha(i)))/sqrt(Nr);
            W_nc_g(:,i) = abs(W_g'*arx);
            W_nc_p(:,i) = abs(W_p'*arx);
        end
        % Calculate norms of columns of both matrices
        W_nc_g_norms = sqrt(sum(W_nc_g.^2));
        W_nc_p_norms = sqrt(sum(W_nc_p.^2));

        %% Non-coherent received data
        y_nc = abs(y(1:M,:)); % Without impairments
        y_nc_norms = sqrt(sum(y_nc.^2)); % Calculate norms of all columns

        %% Non-coherent estimation with imperfect codebooks (without Netwon-Raphson refinement)
        % With gain error
        correlation_gain = ( (W_nc_g./repmat(W_nc_g_norms,M,1))' * (y_nc./repmat(y_nc_norms,M,1)) ); % Size   length(alpha) x training_data_size
        [~, max_corr_indices_gain] = max(correlation_gain); % For each training sample, determine the highest correlation
        AoA_est_gain = alpha(max_corr_indices_gain); % For each training sample, find the best AoA
        % With phase error
        correlation_phase = ( (W_nc_p./repmat(W_nc_p_norms,M,1))' * (y_nc./repmat(y_nc_norms,M,1)) ); % Size   length(alpha) x training_data_size
        [~, max_corr_indices_phase] = max(correlation_phase); % For each training sample, determine the highest correlation
        AoA_est_phase = alpha(max_corr_indices_phase); % For each training sample, find the best AoA

        %% Calculate AoA RMSE
        AoA_RMSE(m_index, error_index, 1) = sqrt(mean(((AoA - AoA_est_gain)/pi*180).^2)); % With gain error
        AoA_RMSE(m_index, error_index, 2) = sqrt(mean(((AoA - AoA_est_phase)/pi*180).^2)); % With phase error

        %% Calculate detection probability (errors smaller than 105/Nr are considered to be correct detection)
        errors_abs_gain = abs((AoA - AoA_est_gain)/pi*180);
        detection_probability(m_index, error_index, 1) = nnz(errors_abs_gain<=(105/Nr))/training_data_size;
        errors_abs_phase = abs((AoA - AoA_est_phase)/pi*180);
        detection_probability(m_index, error_index, 2) = nnz(errors_abs_phase<=(105/Nr))/training_data_size;

    end
    
end


%% Plot results
legend_text = {};
for m_index = 1:length(M_options)
    legend_text{end+1} = sprintf('M = %d', M_options(m_index));
end

% AoA RMSE - gain error
figure
plot(gain_error_options, squeeze(AoA_RMSE(1, :, 1)), 'r-o', 'Linewidth', 2)
hold on
plot(gain_error_options, squeeze(AoA_RMSE(2, :, 1)), 'b-o', 'Linewidth', 2)
plot(gain_error_options, squeeze(AoA_RMSE(3, :, 1)), 'k-o', 'Linewidth', 2)
grid on
xlim([min(gain_error_options) max(gain_error_options)])
xlabel('Gain error std \sigma_A [dB]')
ylabel('AoA RMSE [degree]')
title('AoA RMSE with gain error')
legend(legend_text)
hold off

% Detection probability - gain error
figure
plot(gain_error_options, squeeze(detection_probability(1, :, 1)), 'r-o', 'Linewidth', 2)
hold on
plot(gain_error_options, squeeze(detection_probability(2, :, 1)), 'b-o', 'Linewidth', 2)
plot(gain_error_options, squeeze(detection_probability(3, :, 1)), 'k-o', 'Linewidth', 2)
grid on
xlim([min(gain_error_options) max(gain_error_options)])
xlabel('Gain error std \sigma_A [dB]')
ylabel('Probability')
title('Detection probability with gain error')
legend(legend_text)
hold off

% AoA RMSE - phase error
figure
plot(phase_error_options, squeeze(AoA_RMSE(1, :, 2)), 'r-o', 'Linewidth', 2)
hold on
plot(phase_error_options, squeeze(AoA_RMSE(2, :, 2)), 'b-o', 'Linewidth', 2)
plot(phase_error_options, squeeze(AoA_RMSE(3, :, 2)), 'k-o', 'Linewidth', 2)
grid on
xlim([min(phase_error_options) max(phase_error_options)])
xlabel('Phase error std \sigma_P [deg]')
ylabel('AoA RMSE [degree]')
title('AoA RMSE with phase error')
legend(legend_text)
hold off

% Detection probability - phase error
figure
plot(phase_error_options, squeeze(detection_probability(1, :, 2)), 'r-o', 'Linewidth', 2)
hold on
plot(phase_error_options, squeeze(detection_probability(2, :, 2)), 'b-o', 'Linewidth', 2)
plot(phase_error_options, squeeze(detection_probability(3, :, 2)), 'k-o', 'Linewidth', 2)
grid on
xlim([min(phase_error_options) max(phase_error_options)])
xlabel('Phase error std \sigma_P [degree]')
ylabel('Probability')
title('Detection probability with phase error')
legend(legend_text)
hold off