clear;

%% Directory saves
save_CSV = 1; % flag to save CSV

%save_dir = './data/' % Directory for csv files
save_dir = './data/nn_sim9/'

include_nr_inName = 1; % Include the current Nr in the output file name

%% Basic parameters
Nr = 36; % Rx antenna size %12:12:60
Nt = 1; % Tx antenna size
M = Nr; % number of measurements
d = 0.57; % Array antenna spacing
%SNR_options = 30;%-15:5:30;
SNR_options = [20, 30];
training_data_size = 1e4; % ML training set size
sim_gain_mismatch = 0;
rng(3, 'twister')

%% Generate combining and precoding matrices
sim_exp_bit_match = 0; % 1 or 0
if sim_exp_bit_match
    file_name = ['pn36.csv'];
    steer_vec_rx_phase = csvread(file_name,1,0);
    W = exp(1j*steer_vec_rx_phase(1:M,1:Nr)'/180*pi);
else
    W = ((randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3))/sqrt(2);
end
if Nt > 1
    F = ((randi(2,Nt,M)*2-3)+1j*(randi(2,Nt,M)*2-3))/sqrt(2);
else
    F = 1;
end

%% Generate angles
% AoA = rand(1, training_data_size)*2*pi/4 - pi/4;
angle_region_left = -26.6/180*pi;
angle_region_right = 43.4/180*pi;
angle_region_len = angle_region_right - angle_region_left;
AoA_steps = floor(angle_region_len/(100/180*pi/Nr/2))+1;
% angle_region_bias = (angle_region_right - angle_region_left)/2/180*pi;
AoA = rand(1, training_data_size) * angle_region_len + angle_region_left;
AoD = (Nt~=1)*(rand(1, training_data_size)*2*pi/4 - pi/4);
%alpha = linspace(-pi/4, pi/4, floor(2*Nr*94/105));
alpha = linspace(angle_region_left, angle_region_right, AoA_steps);

%% Generate noiseless data
DFT_quan_bit = 3;
phase_scaler = 2*pi/(2^DFT_quan_bit);
y = zeros(M, training_data_size);
random_phases = rand(1, training_data_size)*pi;
W_DFT = exp(-1j*(0:Nr-1).'*2*pi*d*sin(alpha));
W_DFT_phase = (0:Nr-1).'*2*pi*d*sin(alpha);
W_DFT_phase_quan = round(W_DFT_phase/phase_scaler)*phase_scaler;
W_DFT_quan = exp(-1j*W_DFT_phase_quan);

if sim_gain_mismatch
    W([8,20,32],:) = 2 * W([8,20,32],:);
    W_DFT([8,20,32],:) = 2 * W_DFT([8,20,32],:);
end

for i = 1:training_data_size
    arx = exp(-1j*2*pi*d*(0:Nr-1).'*sin(AoA(i)))/sqrt(Nr);
    atx = exp(-1j*2*pi*d*(0:Nt-1).'*sin(AoD(i)))/sqrt(Nt);
    H = exp(1j*random_phases(i))*arx*atx';
    if Nt > 1
        y(:,i) = diag(W.'*H*F);
    else
        y(:,i) = W.'*H*F;
        
        % DFT beam sounding
        y_DFT(:,i) = W_DFT'*H*F;
        y_DFT_quan(:,i) = W_DFT_quan'*H*F;
    end
end

for snr_index = 1:length(SNR_options)
    
    SNR = SNR_options(snr_index);

    % Add noise
    P_n = 10^(-SNR/10); % Noise power in each measurement, assuming unit signal power
    N = sqrt(P_n/2/Nr)*randn(Nr, training_data_size) + 1j*sqrt(P_n/2/Nr)*randn(Nr, training_data_size);
    y_noisy = y + W'*N; % We consider first M measurements
   
    % Non-coherent received data
    y_nc = abs(y_noisy);
    
    y_DFT_noisy = y_DFT + W_DFT'*N;
    y_DFT_nc = 20*log10(abs(y_DFT_noisy));
    [~,label_DFT_based] = max(y_DFT_nc,[],1);
    
    y_DFT_noisy_quan = y_DFT_quan + W_DFT_quan'*N;
    y_DFT_nc_quan = 20*log10(abs(y_DFT_noisy_quan));
    [~,label_DFT_based_quan] = max(y_DFT_nc_quan,[],1);
    
    if save_CSV
        % save complex measurement
%         if include_nr_inName
%             file_name = [save_dir, 'measurement_raw_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
%         else
%             file_name = [save_dir,'measurement_raw_',num2str(SNR),'dB.csv'];
%         end
%         fprintf('Save file %s \n', file_name)
%         csvwrite(file_name, y_noisy)

        % save RSSI (real non-negative) measurement
        if include_nr_inName
            file_name = [save_dir, 'measurement_RSS_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
        else
            file_name = [save_dir,'measurement_RSS_',num2str(SNR),'dB.csv'];
        end
        fprintf('Save file %s \n', file_name)
        csvwrite(file_name, y_nc)

        % save label
        if include_nr_inName
            file_name = [save_dir, 'label_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
        else
            file_name = [save_dir,'label_',num2str(SNR),'dB.csv'];
        end
        fprintf('Save file %s \n', file_name)
        csvwrite(file_name, AoA)
        fprintf('Done\n')
        
        % save DFT sounding output
        if include_nr_inName
            file_name = [save_dir, 'DFT_output_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
        else
            file_name = [save_dir,'DFT_output_',num2str(SNR),'dB.csv'];
        end
        fprintf('Save file %s \n', file_name)
        csvwrite(file_name, y_DFT_nc)
        fprintf('Done\n')
        
        % save DFT based label
        if include_nr_inName
            file_name = [save_dir, 'DFT_label_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
        else
            file_name = [save_dir,'DFT_label_',num2str(SNR),'dB.csv'];
        end
        fprintf('Save file %s \n', file_name)
        csvwrite(file_name, label_DFT_based)
        fprintf('Done\n')
        
        % save DFT sounding output (quantized DFT beam)
        if include_nr_inName
            file_name = [save_dir, 'DFT_output_quan_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
        else
            file_name = [save_dir,'DFT_output_quan_',num2str(SNR),'dB.csv'];
        end
        fprintf('Save file %s \n', file_name)
        csvwrite(file_name, y_DFT_nc_quan)
        fprintf('Done\n')
        
        % save DFT based label (quantized DFT beam)
        if include_nr_inName
            file_name = [save_dir, 'DFT_label_quan_nr', num2str(Nr), '_', num2str(SNR), 'dB.csv'];
        else
            file_name = [save_dir,'DFT_label_quan',num2str(SNR),'dB.csv'];
        end
        fprintf('Save file %s \n', file_name)
        csvwrite(file_name, label_DFT_based_quan)
        fprintf('Done\n')
        
        
        
        
    end
    
end


