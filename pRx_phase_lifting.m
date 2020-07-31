% Alignment performance evaluation

clear;clc;

%%
% Concantenate all csv files into one variable
x_train = [];
for batch_idx = 1:4
    % Need to manually delete column for 'link_dir', 'tx_beam_label', 'rx_beam_label'
    % 'tx_mask', and 'rx_mask' in this csv files
    
    file_name = ['results_awv0_20-04-15_part',num2str(batch_idx),'_ext.csv'];
    x_temp = csvread(file_name,1,0);
    x_train = [x_train;x_temp];
end

% test data
batch_idx = 5;
file_name = ['results_awv0_20-04-15_part',num2str(batch_idx),'_ext.csv'];
x_test = csvread(file_name,1,0);

x = [x_train;x_test];

% Parameters
target_filter_beam_idx = 35;         % change this to change DFT peak detection target
PN_beam_num = 36;                    % fixed
DFT_beam_num = 64;                   % fixed
axis_offset = 2000;                  % workaround for plotting; don't change this
data_len = size(x,1);                % length of data
train_len = size(x_train,1);         % length of training data
test_len = size(x_test,1);           % length of test data
x_rx_beam = x(1:data_len,2);   % the 4-th column should be rx_beam_idx
x_pRx = x(1:data_len,18);      % the 18-th column should be pRx
cnt = 1;
x_axis_angle = [linspace(-45,45,DFT_beam_num),...
                linspace(axis_offset+1,axis_offset+PN_beam_num,PN_beam_num)];
% figure('Position',[100 100 1200 800])
zz=1;
% legendtext={}
dictionary_train_cnt = zeros(64,1);
dictionary_train = zeros(36,64,100);
pRx_store = ones(64,64,100)*(-101);

% DFT peak detection, filter, and plot
for ii=1:data_len-1
    if x_rx_beam(ii)>x_rx_beam(ii+1)
%         data_to_plot = [x_axis_angle(1+x_rx_beam(cnt:ii))',x_pRx(cnt:ii)];
        
        pRx_this_chunk = x_pRx(cnt:ii);
        rx_beam_this_chunk = x_rx_beam(cnt:ii);

        [~,pRx_max_idx] = max(pRx_this_chunk);
        dict_item_idx = rx_beam_this_chunk(pRx_max_idx);
        PN_beam_idcies = rx_beam_this_chunk(find(rx_beam_this_chunk>63))-63;
        DFT_beam_idcies = rx_beam_this_chunk(find(rx_beam_this_chunk<=63))+1;
        if dict_item_idx<=63
            if logical(PN_beam_idcies)
            dictionary_train_cnt(dict_item_idx+1) = dictionary_train_cnt(dict_item_idx+1) + 1;

            dictionary_train(PN_beam_idcies,dict_item_idx+1,dictionary_train_cnt(dict_item_idx+1)) = ...
                10.^(pRx_this_chunk(find(rx_beam_this_chunk>63))/10);
            pRx_store(DFT_beam_idcies,dict_item_idx+1,dictionary_train_cnt(dict_item_idx+1)) = ...
                pRx_this_chunk(find(rx_beam_this_chunk<=63));
            end
        end

        cnt = ii+1;
    end
    
%     if ii+1 == data_len
%         data_to_plot = [x_axis_angle(1+x_rx_beam(cnt:ii))',x_pRx(cnt:ii)];
%         [~,zzz] = sort(x_axis_angle(1+x_rx_beam(cnt:ii)));
%         cnt = ii;
%     end
end

%% Quick watch of PN beam output
target_filter_beam_idx = 64;
test_data_starts = zeros(64,1);
dict = zeros(36,64);
well_trained = [];
for DFT_beam_idx = 1:64
    if dictionary_train_cnt(DFT_beam_idx)>20
        well_trained = [well_trained,DFT_beam_idx];
        test_data_starts(DFT_beam_idx) = 21;
        ydata = squeeze(dictionary_train(:,DFT_beam_idx,1:20));
        dict(:,DFT_beam_idx) = sum(ydata,2)/20;
    end
end

% Make sure no NaN in this matrix
dict(isnan(dict))=0;

%% Signal processing parameters

M = 64;                 % number of measurement
MCS = M/2;               % effective number of measurement in complex CS problem
Nr = 36;
pattern_watch = 0; % 1 is PN, 0 is DFT

if pattern_watch
    % load PN from csv
    M = 30;                 % number of measurement
    MCS = M/2;               % effective number of measurement in complex CS problem
    file_name = ['pn36.csv'];
    steer_vec_rx_phase = csvread(file_name,1,0);
    steer_vec_rx = exp(1j*steer_vec_rx_phase(1:M,:)'/180*pi);
else
    file_name = ['dft64.csv'];
    steer_vec_rx_phase = csvread(file_name,1,0);
    steer_vec_rx = exp(1j*steer_vec_rx_phase(1:M,:)'/180*pi);
end

%% dictionary generation
cand_num_r = 64;
dict_pRx_norm = zeros(cand_num_r,1);

cand_y = zeros(M,cand_num_r);
cand_angle_r = linspace(-pi*45/180,pi*45/180,cand_num_r);
AOAstep = cand_angle_r(2) - cand_angle_r(1);

cand_ARV_r = exp(-1j*(0:Nr-1)'*0.56*2*pi*sin(cand_angle_r)); % seems d = 0.56lambda works best
%% quick watch on actual PN beam pattern and simulated ones
B_index_target = 21;
y_train = dict(B_index_target,14:end).';
new_dict = cand_ARV_r(:,14:end);
A_train = new_dict.';
pattern = abs(A_train * steer_vec_rx(:,B_index_target)).^2;
figure
target_filter_beam_idx = B_index_target;

if pattern_watch
    all_pattern = squeeze(dictionary_train(target_filter_beam_idx,:,1:20));  
else
    all_pattern = 10.^(squeeze(pRx_store(target_filter_beam_idx,:,1:20))/10);
end

% linear scale with mean and error bar with std. dev.
% the issues is details of beam patterns are missing 
%     all_pattern = squeeze(dictionary_train(target_filter_beam_idx,:,1:20));
%     mean_pattern = sum(all_pattern,2)/20;
%     plot(cand_angle_r(14:end)/pi*180,1e8*(mean_pattern(14:end)),'k-','linewidth',3)
%     hold on
%     temp = sqrt(var(all_pattern'));
%     er = errorbar(cand_angle_r(14:end)/pi*180,...
%                   1e8*(mean_pattern(14:end)),...
%                   1e8*temp(14:end));   
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';  


% dB scale with mean and error bar for value with certain percentile
mean_pattern = sum(all_pattern,2)/20;
plot(cand_angle_r(14:end)/pi*180,10*log10(mean_pattern(14:end)),'k-','linewidth',3)
hold on
temp = sort(all_pattern,2);
er = errorbar(cand_angle_r(14:end)/pi*180,...
              10*log10(mean_pattern(14:end)),...
              10*log10(mean_pattern(14:end)) - 10*log10(temp(14:end,2)),...
              -10*log10(mean_pattern(14:end)) + 10*log10(temp(14:end,18)));   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  


hold on
grid on
xlabel('Angle [deg]')
ylabel('Received Power [dB]')
xlim([-45,45])
% ylim([-100,-55])
% title('Measured Pattern w/ PN Codebook')

MMSE_scaling = pinv(pattern)*mean_pattern(14:end);
plot(cand_angle_r(14:end)/pi*180,10*log10(MMSE_scaling*pattern),'b-','linewidth',3)
hold on
% title('Simulated Pattern w/ PN Codebook')
% grid on
% ylabel('Normalized pRx [dB]')
xlim([-30,50])
% ylim([-2,82])
ylim([-100,-55])
set(gca,'FontSize',14)
legend('Measured (Mean)','Measured (Std. Dev.)','Model')
%%
mainlobe_area = B_index_target-13-1:B_index_target-13+1;
sidelobe_area = [1:B_index_target-13-2,B_index_target-13+2:51];
alllobe_area = 1:51;

pattern1 = MMSE_scaling*pattern;
pattern2 = mean_pattern(14:end);
lobe_area = mainlobe_area;
distortion_metric = 20*log10(norm(pattern1(lobe_area) - pattern2(lobe_area))/norm(pattern2(lobe_area)));
distortion_metric
% print(distortion_metric);
%% Refine knowledge of APR using training data
        
[Ut, Sigmat, Vt] = svd(steer_vec_rx.');

APR_raw = Ut(:,1:MCS)*sqrt(Sigmat(1:MCS,1:MCS));
ACS = sqrt(Sigmat(1:MCS,1:MCS))*Vt(:,1:MCS)';

new_dict = cand_ARV_r(:,14:end); 

dict_pRx = abs(steer_vec_rx.'*cand_ARV_r);

Measure_mat_new = ACS*cand_ARV_r;
for cc=1:cand_num_r
    if cc>=14
        Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
    else
        Measure_mat_new_norm(cc) = 1e6; % force system not pick them
    end
end

% Estimate dictionary APR using PhaseLift
A_train = new_dict.';
for mm=1:M
    % y^H = PN^H(:,mm) * A_rx^H
    % y = A_rx * ACS^H * APR^H(:,mm)
    % where PN = APR * ACS

    y_train = dict(mm,14:end).' * 1e8; % PN pattern from measurement; Scaling to around 1
    AthetACS = A_train * ACS.';
    cvx_begin sdp quiet
        variable ZZ(MCS,MCS) hermitian semidefinite 
        minimize 1*trace(ZZ)+ norm(diag(AthetACS * ZZ * AthetACS') - y_train, 2)
    cvx_end
    [Umat,Sigma,Vmat] = svd(ZZ);
    APR_hat(:,mm) = Vmat(:,1);

end


%% Process actual data
confusion_mtx = zeros(64,64);
confusion_mtx_prob = zeros(64,64);
gain_loss = ones(64,80)*(-99);

for DFT_beam_idx = 1:64
    if logical(find(well_trained==DFT_beam_idx)>0)
        ydata = squeeze(dictionary_train(:,DFT_beam_idx,test_data_starts(DFT_beam_idx):end));
        gain_ydata = squeeze(pRx_store(:,DFT_beam_idx,test_data_starts(DFT_beam_idx):end));
        
        
        for ii=1:dictionary_train_cnt(DFT_beam_idx) - test_data_starts(DFT_beam_idx)+1
            score = zeros(64,1);

            % Normalize such that it is around 1
            y_data_norm = ydata(1:M,ii)*1e8;
            
            % Make sure no NaN
            y_data_norm(isnan(y_data_norm))=0;
            
            % Use data enhanced APR
            APR = APR_hat.';
            
            % Use raw APR (directly from SVD of W matrix)
%             APR = APR_raw;

            % PhaseLift Main Solver
            clc;fprintf('Evaluate data %d in True AoA Index %d\n',ii,DFT_beam_idx);
            cvx_begin sdp quiet
                variable Z(MCS,MCS) hermitian semidefinite 
                minimize trace(Z) + norm(diag(APR * Z * APR') - y_data_norm, 2)
            cvx_end

            % Getting first eigenvector as solution of PhaseRetrieval problem
            [Umat,Sigma,Vmat] = svd(Z);      
        
            for cand_idx = 1:64
                if logical(find(well_trained==cand_idx)>0)
                    score(cand_idx) = (Measure_mat_new(:,cand_idx)'*Vmat(:,1))./Measure_mat_new_norm(cand_idx)';
                else
                    score(cand_idx) = 0;
                end
            end
            % update beam alignment matrix
            [~,dict_item_best] = max(score);
            confusion_mtx(DFT_beam_idx, dict_item_best) = confusion_mtx(DFT_beam_idx, dict_item_best)+1;

            % update gain matrix
            gain_loss(DFT_beam_idx,ii) = max(gain_ydata(:,ii)) - gain_ydata(dict_item_best,ii);
        end
    end
end


for DFT_beam_idx = 1:64
    if sum(confusion_mtx(DFT_beam_idx,:))>0
        confusion_mtx_prob(DFT_beam_idx,:) = confusion_mtx(DFT_beam_idx,:)/sum(confusion_mtx(DFT_beam_idx,:));
    end
end

% Confusion matrix
figure
X = 1:64;
Y = 1:64;
[xplot,yplot] = meshgrid(X,Y);
figure('Position',[100 100 800 800])
heatmap(X,Y,confusion_mtx_prob)
colorbar

det_prob = sum(diag(confusion_mtx))/sum(sum(confusion_mtx));
%%
for PN_size_idx = 1:1
    temp = reshape(gain_loss,80*64,1);
    gain_loss_this_PN_num = temp(find(temp>-0.01));
    temp2 = sort(gain_loss_this_PN_num,'ascend');
    u_idx = floor(length(temp2)*0.9);
    l_idx = floor(length(temp2)*0.1);
    
    bar_mean = mean(gain_loss_this_PN_num);
    bar_low = bar_mean - temp2(l_idx);
    bar_high = -bar_mean + temp2(u_idx);
end

% xdata = 1:length(PN_size_range)*5;
lightgreen = [0.75, 0.75, 0];
green = 	[0.4660, 0.6740, 0.1880];

barWidth = 0.2;
figure('position',[100,100,600,500])
b = bar(1,[bar_mean;bar_high],'stacked','BarWidth', barWidth) 
b(1).FaceColor = lightgreen;
b(2).FaceColor = green;
hold on
% bar(xdata(1:5:end),[bar_mean-bar_mean;bar_high-bar_high],'stacked') 
% hold on
% bar(xdata(5:5:end),[bar_mean-bar_mean;bar_high-bar_high],'stacked') 
% hold on
% grid on
% xticks(xdata(3:5:end))
% myxtick = [sprintf('% 2d',PN_size_range(1))];
% for ii=1:length(PN_size_range)-1
%     if PN_size_range(ii+1)<10
%         myxtick = [myxtick;sprintf('% 2d',PN_size_range(ii+1))];
%     else
%         myxtick = [myxtick;sprintf('%2d',PN_size_range(ii+1))];
%     end
% end
% set(gca,'xticklabel',myxtick)
% er = errorbar(PN_size_range,...
%               bar_mean,...
%               bar_low,...
%               bar_high);   
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% ylim([0,20])
grid on
set(gca,'FontSize',14)
ylabel('Beamforming Gain Loss [dB]')
xlabel('Number of Measurements (M)')
ylim([0,15]) 
% xlim([1,max(xdata)-1])
legend('PhaseLift w/ Dict. Ref., (50% percentile)',...
       'PhaseLift w/ Dict. Ref., (90% percentile)')

% %%
% 
% 
% for runindex=1:runtimes
%     clc;fprintf('Ite %d out of %d\n',runindex,runtimes);
%            
%     APR = APR_hat.';
% %         y = (APR * ACS * new_dict(:,32)) + noise_normal*sqrt(noise_pow);
%     y = (steer_vec_rx.' * new_dict(:,32)) + noise_normal*sqrt(noise_pow);
% 
%     d_vec = (abs(y)).^2;
% 
%     cvx_begin sdp quiet
%         variable Z(MCS,MCS) hermitian semidefinite 
%         minimize trace(Z) + norm(diag(APR * Z * APR') - d_vec,2)
%     cvx_end
% 
%     z_new = (ACS*new_dict(:,32));
%     Z_new = z_new*z_new';
%     norm(diag(APR * Z_new * APR') - d_vec);
% 
%     [Umat,Sigma,Vmat] = svd(Z);
%     alpha = pinv(Vmat(:,1))*z_new;
%     norm(alpha*Vmat(:,1)-z_new)/norm(z_new);
% 
%     % PhaseLift Solution
%     cand_score_PLCS = (Measure_mat_new'*Vmat(:,1))./Measure_mat_new_norm';
%     [~,bestindex_PLCS(runindex)] = max(abs(cand_score_PLCS));
%     bestAOA_PLCS(runindex) = (bestindex_PLCS(runindex)-1)*AOAstep-45*pi/180;
%     AOA_error_PLCS(runindex) = abs(bestAOA_PLCS(runindex) - 0);
% 
% 
% end
% %%
% AOAalign_PLCS_mean = sum((AOA_error_PLCS/pi*180)<(105/Nr),1)/runtimes;
% 
% figure
% plot(SNR_range,AOAalign_PLCS_mean,'-','linewidth',2);
% grid on
% xlabel('SNR (dB)')
% ylabel('Alignment Prob')
% legend('PhaseLift')

