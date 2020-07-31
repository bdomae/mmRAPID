clear;
run('data_for_compressive_sensing')


dict_type = 0; % 1 is model predicted; 0 is estimated from training data
%% Basic parameters
% alpha = linspace(-pi/4, pi/4, floor(2*Nr*90/105));
M_options = 4:1:12;
% M0 = max(M_options);

%% Generate non-coherent codebook
W_nc_original = zeros(M, length(alpha));
for i = 1:length(alpha)
    arx = exp(-1j*(0:Nr-1)'*pi*sin(alpha(i)))/sqrt(Nr);
    W_nc_original(:,i) = abs(W'*arx);
end
data_size = length(AoA);

%% AoA estimation
y_original = y; % Noiseless original
AoA_RMSE = zeros(length(M_options), 1);
detection_probability = zeros(length(M_options), 1);
legend_text = {};
for snr_index = 1
    SNR = SNR_options(snr_index);
%         SNR = SNR_options(snr_index);

    % Read CSV of raw data
    filename = [save_dir(1:end),'measurement_RSS_',num2str(SNR),'dB.csv'];
    y_nc_all_data = csvread(filename);
    M_max = size(y_nc_all_data,1);
    
    % Read CSV of label (AoA)
    filename = [save_dir(1:end),'label_',num2str(SNR),'dB.csv'];
    AoA = csvread(filename);
    
    % Dictionary estimation using DFT labels
    train_cnt = zeros(length(alpha),1);
    dict_est_sum = zeros(M_max,length(alpha));
    traing_size = 7000;
    for dd=1:traing_size
        best_DFT_bin = label_DFT_based(dd);
        train_cnt(best_DFT_bin) = train_cnt(best_DFT_bin) + 1;
        dict_est_sum(:,best_DFT_bin) = dict_est_sum(:,best_DFT_bin) + y_nc_all_data(:,dd);
    end
    dict_est = dict_est_sum./repmat(train_cnt.',M_max,1);
    
    for m_index = 1:length(M_options)
        M = M_options(m_index);
        legend_text{end+1} = sprintf('M = %d meas.', M);
        
        if dict_type == 1
            W_nc = W_nc_original(1:M,:); % Select subset of non-coherent combiners
        else
            W_nc = dict_est(1:M,:);
        end
        W_nc_norms = sqrt(sum(W_nc.^2)); % Calculate norms of all columns   
        
        y_nc = y_nc_all_data(1:M,traing_size+1:end); 
        label_DFT_based_test = label_DFT_based(traing_size+1:end);
        test_data_size = data_size - traing_size;
        y_DFT_nc_test = y_DFT_nc(:,traing_size+1:end);
       
        %% Non-coherent received data
%         y_nc = abs(y(1:M,:));
        y_nc_norms = sqrt(sum(y_nc.^2)); % Calculate norms of all columns 

        %% Non-coherent estimation (without Netwon-Raphson refinement)
        correlation = ( (W_nc./repmat(W_nc_norms,M,1))' * (y_nc./repmat(y_nc_norms,M,1)) ); % Size   length(alpha) x training_data_size
        [~, max_corr_indices] = max(correlation); % For each training sample, determine the highest correlation
        AoA_est = alpha(max_corr_indices); % For each training sample, find the best AoA

        %% Calculate AoA RMSE
        AoA_RMSE(m_index, snr_index) = sqrt(mean(((AoA(traing_size+1:end) - AoA_est)/pi*180).^2));
%         gain_loss = 10*log10(Nr^2)*ones(1,data_size) - 10*log10(abs(sum(exp(1j*pi*(0:Nr-1).'*sin(AoA - AoA_est)),1)).^2);
        for dd = 1:test_data_size
            gain_loss(dd) = y_DFT_nc_test(label_DFT_based_test(dd),dd) - y_DFT_nc_test(max_corr_indices(dd),dd);
        end
        temp = sort(gain_loss,'ascend');
        gain_loss_mean(m_index, snr_index) = temp(floor(test_data_size*0.5));
        gain_loss_high(m_index, snr_index) = temp(floor(test_data_size*0.9));
        
        %% Calculate detection probability (errors smaller than 105/Nr are considered to be correct detection)
        errors_abs = abs((AoA(traing_size+1:end) - AoA_est)/pi*180);
        detection_probability(m_index, snr_index) = nnz(errors_abs<=(105/Nr))/test_data_size;
        
    end
    
end

%% Plot results

figure
plot(M_options,gain_loss_high,'-o','linewidth',2)
hold on
plot(M_options,gain_loss_mean,'-o','linewidth',2)
grid on
yticks([1,2,3,4,5,10,20])
xlabel('Number of Measurement (M)')
ylabel('Gain Loss [dB]')
legend('90 Percentile','50 Percentile')
%%
% AoA RMSE
figure
plot(M_options, AoA_RMSE, 'r-o', 'Linewidth', 2)
hold on

grid on
set(gca,'FontSize',14)
% xlim([min(SNR_options) max(SNR_options)])
% xticks(SNR_options)
xlabel('Number of Measurement (M)')
ylabel('AoA RMSE [degree]')
title('RMSE vs SNR')
legend(legend_text)
hold off

% Detection probability
figure
plot(M_options, detection_probability, 'r-o', 'Linewidth', 2)
grid on
set(gca,'FontSize',14)
% xticks(SNR_options)
% xlim([min(SNR_options) max(SNR_options)])
xlabel('Number of Measurement (M)')
ylabel('Alignment probability')
title('Alignment Probability vs SNR')
legend(legend_text)
hold off
%% save results to CSV
% file_name = [save_dir,'RSSI_MP_RMSE.csv'];
% fprintf('Save file %s \n', file_name)
% csvwrite(file_name(2:end), AoA_RMSE)
% 
% file_name = [save_dir,'RSSI_MP_alignment_rate.csv'];
% fprintf('Save file %s \n', file_name)
% csvwrite(file_name(2:end), detection_probability)
% fprintf('Done\n')
%%
lightgreen = [0.0431, 0.4000, 0.1373];
blue = 	[0, 0.4470, 0.7410];
x=12:12:72;
x_theo=8:2:108;
figure;
loglog(x_theo,2*x_theo,'k-','linewidth',2);
hold on
grid on;
loglog(x,[6,7,7,7,8,8],'s','linewidth',2,'markersize',12,'Color',blue);
hold on;
loglog(x_theo,1.4*log2(x_theo),'-','linewidth',2,'Color',blue);
hold on;
loglog([36],[20],'^','linewidth',2,'markersize',12,'Color',blue);
hold on;
loglog(x,[4,5,5,7,7,7],'ro','linewidth',2,'markersize',10,'Color',lightgreen);
hold on;
loglog(x_theo,1.1*log2(x_theo),'r--','linewidth',2,'Color',lightgreen);
hold on;
loglog([36],[5],'rx','linewidth',2,'markersize',10,'Color',lightgreen);
xticks([12,36,72])
xlim([8,108])
ylim([3,128])
yticks([4,8,16,32,64])
legend({'Exhaustive Search',...
       'RSS-MP (Sim.)',...
       'RSS-MP (Sim. Fitting)',...
       'RSS-MP (Exp.)',...
       'NN (Sim.)',...
       'NN (Sim. Fitting)',...
       'NN (Exp.)'}, 'Location','north','NumColumns',2);
set(gca,'FontSize',13)
xlabel('Antenna Number (N_R)');
ylabel('Required Measurement (M)')
%%
% SNR = 10;
% P_n = 10^(-SNR/10);
% W_nc_original = zeros(1, length(alpha));
% MCtimes = 1e3;
% for MCidx = 1:MCtimes
%     N = sqrt(P_n/2/Nr)*randn(Nr, 1) + 1j*sqrt(P_n/2/Nr)*randn(Nr, 1);
% 
%     for i = 1:length(alpha)
%         arx = exp(-1j*(0:Nr-1)'*pi*sin(alpha(i)))/sqrt(Nr);
%         W_nc_original(MCidx,i) = 20*log10(abs(W(:,1)'*arx + W(:,1)'*N));
%     end
% end
% for rr=1:length(alpha)
%     temp = sort(W_nc_original(:,rr),'ascend');
%     pattern_mean(rr) = temp(MCtimes*0.5);
%     pattern_low(rr) = temp(MCtimes*0.5) - temp(MCtimes*0.1);
%     pattern_high(rr) = -temp(MCtimes*0.5) + temp(MCtimes*0.9);
% end
cnt = 1;
for dd=1:data_size
    if label_DFT_based(dd) == 10
        pattern(:,cnt) = 20*log10(y_nc_all_data(1:M,dd));
        cnt = cnt + 1;
    end
end
pattern_mean = zeros(M,1);
pattern_low = zeros(M,1);
pattern_high = zeros(M,1);
for mm=1:M
    temp = sort(pattern(mm,:),'ascend');
    pattern_mean(mm) = temp(floor(cnt*0.5));
    pattern_low(mm) = temp(floor(cnt*0.5)) - temp(floor(cnt*0.1));
    pattern_high(mm) = -temp(floor(cnt*0.5)) + temp(floor(cnt*0.9));
end

figure
plot(1:M,pattern_mean,'linewidth',3)
hold on
e = errorbar(1:M,...
         pattern_mean,...
         pattern_low,...
         pattern_high);
e.Color = 'black'
grid on
title('Pattern (Sim) with "30dB" SNR')