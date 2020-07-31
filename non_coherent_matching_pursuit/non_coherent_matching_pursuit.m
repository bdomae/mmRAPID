clear;
run('data_for_compressive_sensing')

%% Basic parameters
alpha = linspace(-pi/3, pi/3, 1024);
M_options = 5:5:15;

%% Generate non-coherent codebook
W_nc_original = zeros(M, length(alpha));
for i = 1:length(alpha)
    arx = exp(-1j*(0:Nr-1)'*pi*sin(alpha(i)))/sqrt(Nr);
    W_nc_original(:,i) = abs(W'*arx);
end


%% AoA estimation
y_original = y; % Noiseless original
AoA_RMSE = zeros(length(M_options), length(SNR_options));
detection_probability = zeros(length(M_options), length(SNR_options));
legend_text = {};
for m_index = 1:length(M_options)
    M = M_options(m_index);
    legend_text{end+1} = sprintf('M = %d meas.', M);
    
    W_nc = W_nc_original(1:M,:); % Select subset of non-coherent combiners
    W_nc_norms = sqrt(sum(W_nc.^2)); % Calculate norms of all columns
    
    for snr_index = 1:length(SNR_options)
        
        SNR = SNR_options(snr_index);
        
        % Read CSV of raw data
        filename = [save_dir(2:end),'measurement_raw_',num2str(SNR),'dB.csv'];
        y = csvread(filename);
        
        % Read CSV of label (AoA)
        filename = [save_dir(2:end),'label_',num2str(SNR),'dB.csv'];
        AoA = csvread(filename);

        %% Non-coherent received data
        y_nc = abs(y(1:M,:));
        y_nc_norms = sqrt(sum(y_nc.^2)); % Calculate norms of all columns 

        %% Non-coherent estimation (without Netwon-Raphson refinement)
        correlation = ( (W_nc./repmat(W_nc_norms,M,1))' * (y_nc./repmat(y_nc_norms,M,1)) ); % Size   length(alpha) x training_data_size
        [~, max_corr_indices] = max(correlation); % For each training sample, determine the highest correlation
        AoA_est = alpha(max_corr_indices); % For each training sample, find the best AoA

        %% Calculate AoA RMSE
        AoA_RMSE(m_index, snr_index) = sqrt(mean(((AoA - AoA_est)/pi*180).^2));
        
        %% Calculate detection probability (errors smaller than 105/Nr are considered to be correct detection)
        errors_abs = abs((AoA - AoA_est)/pi*180);
        detection_probability(m_index, snr_index) = nnz(errors_abs<=(105/Nr))/training_data_size;
        
    end
    
end

%% Plot results

% AoA RMSE
figure
plot(SNR_options, AoA_RMSE(1,:), 'r-o', 'Linewidth', 2)
hold on
plot(SNR_options, AoA_RMSE(2,:), 'b-o', 'Linewidth', 2)
plot(SNR_options, AoA_RMSE(3,:), 'k-o', 'Linewidth', 2)
grid on
set(gca,'FontSize',14)
xlim([min(SNR_options) max(SNR_options)])
xticks(SNR_options)
xlabel('SNR [dB]')
ylabel('AoA RMSE [degree]')
title('RMSE vs SNR')
legend(legend_text)
hold off

% Detection probability
figure
plot(SNR_options, detection_probability(1,:), 'r-o', 'Linewidth', 2)
hold on
plot(SNR_options, detection_probability(2,:), 'b-o', 'Linewidth', 2)
plot(SNR_options, detection_probability(3,:), 'k-o', 'Linewidth', 2)
grid on
set(gca,'FontSize',14)
xticks(SNR_options)
xlim([min(SNR_options) max(SNR_options)])
xlabel('SNR [dB]')
ylabel('Alignment probability')
title('Alignment Probability vs SNR')
legend(legend_text)
hold off
%% save results to CSV
file_name = [save_dir,'RSSI_MP_RMSE.csv'];
fprintf('Save file %s \n', file_name)
csvwrite(file_name(2:end), AoA_RMSE)

file_name = [save_dir,'RSSI_MP_alignment_rate.csv'];
fprintf('Save file %s \n', file_name)
csvwrite(file_name(2:end), detection_probability)
fprintf('Done\n')