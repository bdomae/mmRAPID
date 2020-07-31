clear;clc;

% Concantenate all csv files into one variable
x_train = [];
for batch_idx = 1:15
    % Need to manually delete column for 'link_dir', 'tx_beam_label', 'rx_beam_label'
    % 'tx_mask', and 'rx_mask' in this csv files
    % Apr15, tx_gain_idx=15; Apr16, tx_gain_idx=21; Apr20, tx_gain_idx=12 
    if batch_idx <= 5
        data_num = 16;
    elseif batch_idx <= 10
        data_num = 15;
    else
        data_num = 20;
    end
%     file_name = ['results_awv0_20-04-15_part',num2str(batch_idx),'_ext.csv'];
    file_name = ['results_awv0_20-04-',num2str(data_num),'_part',num2str(1+mod(batch_idx-1,5)),'_ext.csv'];
    x_temp = csvread(file_name,1,0);
    x_train = [x_train;x_temp];
end

% test data (uses other way to divide training and test)
% batch_idx = 5;
% % file_name = ['results_awv0_20-04-15_part',num2str(batch_idx),'_ext.csv'];
% file_name = ['results_awv0_20-04-16_part',num2str(batch_idx),'_ext.csv'];
% x_test = csvread(file_name,1,0);

x = x_train;

% Parameters
target_filter_beam_idx = 35;         % change this to change DFT peak detection target
PN_beam_num = 36;                    % fixed
DFT_beam_num = 64;                   % fixed
data_len = size(x,1);                % length of data
train_len = 60;                      % requried training in each DFT bin
x_rx_beam = x(1:data_len,2);   % the 4-th column should be rx_beam_idx
x_pRx = x(1:data_len,18);      % the 18-th column should be pRx


% legendtext={}
dictionary_train_cnt = zeros(64,1);
dictionary_train = zeros(36,64,100);
pRx_store = ones(64,64,200)*(-100);

% DFT peak detection, filter, and plot
cnt = 1;
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
    
end
%% Dictionary directly from model
file_name = ['pn36.csv'];
steer_vec_rx_phase = csvread(file_name,1,0);
M = 36;
Nr = 36;
steer_vec_rx = exp(1j*steer_vec_rx_phase(1:M,:)'/180*pi);

% TG possibly has some gain mismatch in references port (subtle difference)
% steer_vec_rx([8,20,32],:) = 1.5*steer_vec_rx([8,20,32],:);

cand_num_r = 64;
dict_pRx_norm = zeros(cand_num_r,1);

cand_y = zeros(M,cand_num_r);
cand_angle_r = linspace(-pi*44.8/180,pi*43.4/180,cand_num_r);
AOAstep = cand_angle_r(2) - cand_angle_r(1);

cand_ARV_r = exp(-1j*(0:Nr-1)'*0.57*2*pi*sin(cand_angle_r));

new_dict = cand_ARV_r;
A_train = new_dict.';
pattern = (abs(A_train * steer_vec_rx).^2)';


%% Quick watch of DFT beam output
target_filter_beam_idx = 24;
figure
plot(1:64,squeeze(pRx_store(:,target_filter_beam_idx,:)))
hold on
grid on
xlabel('DFT Beam Index')
ylabel('pRx [non-log scale]')
%% Quick watch of PN beam output
target_filter_beam_idx = 64;
test_data_starts = zeros(64,1);
dict = zeros(36,64);
well_trained = [];

% This manualy rule out last DFT beam in training and testing!!!
dictionary_train_cnt(64) = 0;

for DFT_beam_idx = 1:64
    if dictionary_train_cnt(DFT_beam_idx)>train_len
        well_trained = [well_trained,DFT_beam_idx];
        test_data_starts(DFT_beam_idx) = train_len+1;
        ydata = squeeze(dictionary_train(:,DFT_beam_idx,1:train_len));
        dict(:,DFT_beam_idx) = sum(ydata,2)/train_len;
    end
end
PN_pattern = zeros(36,100);
PN_pattern_norm = zeros(36,100);
PN_pattern = squeeze(dictionary_train(:,target_filter_beam_idx,:));
for ii=1:100
    if norm(PN_pattern(:,ii))>0
        PN_pattern_norm(:,ii) = PN_pattern(:,ii)./norm(PN_pattern(:,ii));
    end
end
figure
plot(1:36,PN_pattern_norm)
hold on
% plot(1:36,dict(:,target_filter_beam_idx),'ko-','linewidth',4)
% grid on
xlabel('PN Beam Index')
ylabel('pRx [non-log scale]')
%% quick look at dictionary (from model vs from training data)
angle_of_interest = 35;
figure
subplot(211)
plot(1:36,dict(1:M,angle_of_interest))
grid on
subplot(212)
plot(1:36,pattern(:,angle_of_interest))
grid on
%% Testing w/o Non-dictionary enhanced 
PN_size_range = [4,5,6,7,10,12,15,20];
gain_loss = ones(64,200-train_len,length(PN_size_range))*(-99);
dict_type = 1; % 0 is purely based on model; 1 is based on training data;

for PN_size_idx = 1:length(PN_size_range)
    confusion_mtx = zeros(64,64);
    confusion_mtx_prob = zeros(64,64);
    PN_size = PN_size_range(PN_size_idx);

    for DFT_beam_idx = 1:64
        if logical(find(well_trained==DFT_beam_idx)>0)
            ydata = squeeze(dictionary_train(:,DFT_beam_idx,test_data_starts(DFT_beam_idx):end));
            gain_ydata = squeeze(pRx_store(:,DFT_beam_idx,test_data_starts(DFT_beam_idx):end));

            for ii=1:dictionary_train_cnt(DFT_beam_idx) - test_data_starts(DFT_beam_idx)+1
                score = zeros(64,1);
                for cand_idx = 1:64
                    if logical(find(well_trained==cand_idx)>0)
                        if dict_type
                            score(cand_idx) = abs(ydata(1:PN_size,ii)'*dict(1:PN_size,cand_idx)/norm(dict(1:PN_size,cand_idx)));
                        else
                            score(cand_idx) = abs(ydata(1:PN_size,ii)'*pattern(1:PN_size,cand_idx)/norm(pattern(1:PN_size,cand_idx)));
                        end
                    else
                        score(cand_idx) = 0;
                    end
                end
                % update beam alignment matrix
                [~,dict_item_best] = max(score);
                confusion_mtx(DFT_beam_idx, dict_item_best) = confusion_mtx(DFT_beam_idx, dict_item_best)+1;

                % update gain matrix
                gain_loss(DFT_beam_idx,ii,PN_size_idx) = max(gain_ydata(:,ii)) - gain_ydata(dict_item_best,ii);
            end
        end
    end


    for DFT_beam_idx = 1:64
        if sum(confusion_mtx(DFT_beam_idx,:))>0
            confusion_mtx_prob(DFT_beam_idx,:) = confusion_mtx(DFT_beam_idx,:)/sum(confusion_mtx(DFT_beam_idx,:));
        end
    end

    % X = 1:64;
    % Y = 1:64;
    % [xplot,yplot] = meshgrid(X,Y);
    % figure('Position',[100 100 800 800])
    % heatmap(X,Y,confusion_mtx_prob)
    % colorbar

    det_prob(PN_size_idx) = sum(diag(confusion_mtx))/sum(sum(confusion_mtx));
end

%% Testing w/ Non-dictionary enhanced
PN_size_range = [4,5,6,7,10,12,15,20];
gain_loss2 = ones(64,200-train_len,length(PN_size_range))*(-99);
dict_type = 0; % 0 is purely based on model; 1 is based on training data;

for PN_size_idx = 1:length(PN_size_range)
    confusion_mtx = zeros(64,64);
    confusion_mtx_prob = zeros(64,64);
    PN_size = PN_size_range(PN_size_idx);

    for DFT_beam_idx = 1:64
        if logical(find(well_trained==DFT_beam_idx)>0)
            ydata = squeeze(dictionary_train(:,DFT_beam_idx,test_data_starts(DFT_beam_idx):end));
            gain_ydata = squeeze(pRx_store(:,DFT_beam_idx,test_data_starts(DFT_beam_idx):end));

            for ii=1:dictionary_train_cnt(DFT_beam_idx) - test_data_starts(DFT_beam_idx)+1
                score = zeros(64,1);
                for cand_idx = 1:64
                    if logical(find(well_trained==cand_idx)>0)
                        if dict_type
                            score(cand_idx) = abs(ydata(1:PN_size,ii)'*dict(1:PN_size,cand_idx)/norm(dict(1:PN_size,cand_idx)));
                        else
                            score(cand_idx) = abs(ydata(1:PN_size,ii)'*pattern(1:PN_size,cand_idx)/norm(pattern(1:PN_size,cand_idx)));
                        end
                    else
                        score(cand_idx) = 0;
                    end
                end
                % update beam alignment matrix
                [~,dict_item_best] = max(score);
                confusion_mtx(DFT_beam_idx, dict_item_best) = confusion_mtx(DFT_beam_idx, dict_item_best)+1;

                % update gain matrix
                gain_loss2(DFT_beam_idx,ii,PN_size_idx) = max(gain_ydata(:,ii)) - gain_ydata(dict_item_best,ii);
            end
        end
    end


    for DFT_beam_idx = 1:64
        if sum(confusion_mtx(DFT_beam_idx,:))>0
            confusion_mtx_prob(DFT_beam_idx,:) = confusion_mtx(DFT_beam_idx,:)/sum(confusion_mtx(DFT_beam_idx,:));
        end
    end
    
%     % Confusion matrix with this M number
%     X = 1:64;
%     Y = 1:64;
%     [xplot,yplot] = meshgrid(X,Y);
%     figure('Position',[100 100 800 800])
%     heatmap(X,Y,confusion_mtx_prob)
%     colorbar

    det_prob(PN_size_idx) = sum(diag(confusion_mtx))/sum(sum(confusion_mtx));
end
%%
high_percentile = 0.9;
for PN_size_idx = 1:length(PN_size_range)
    temp = reshape(squeeze(gain_loss(:,:,PN_size_idx)),(200-train_len)*64,1);
    gain_loss_this_PN_num = temp(find(temp>-0.01));
    temp2 = sort(gain_loss_this_PN_num,'ascend');
    u_idx = floor(length(temp2)*high_percentile);
    l_idx = floor(length(temp2)*(1-high_percentile));
    
    bar_mean(PN_size_idx) = mean(gain_loss_this_PN_num);
    bar_low(PN_size_idx) = temp2(l_idx);
    bar_high(PN_size_idx) = temp2(u_idx);
end

%%
for PN_size_idx = 1:length(PN_size_range)
    temp = reshape(squeeze(gain_loss2(:,:,PN_size_idx)),(200-train_len)*64,1);
    gain_loss_this_PN_num = temp(find(temp>-0.01));
    temp2 = sort(gain_loss_this_PN_num,'ascend');
    u_idx = floor(length(temp2)*high_percentile);
    l_idx = floor(length(temp2)*(1-high_percentile));
    
    bar_mean2(PN_size_idx) = mean(gain_loss_this_PN_num);
    bar_low2(PN_size_idx) = temp2(l_idx);
    bar_high2(PN_size_idx) = temp2(u_idx);
end
%%
figure
plot(PN_size_range,1-det_prob)
xlabel('Number of Beacon')
ylabel('Mis-Alignment Probability')
grid on

%%
% NN data (from TF analysis) PN_size_range = [4,5,6,7,10,12,15];
bar_mean_nn = [0, 0, 0, 0, 0, 0, 0,0];
bar_high_nn = [19.92, 1.15, 0.21, 0, 0, 0, 0,0];

% Bar graph creation
xdata = 1:length(PN_size_range)*5;
lightblue = [0.0 0.8 1];
blue = 	[0, 0.4470, 0.7410];
lightred = [0.9, 0.5, 0.0];
red = [0.8500, 0.3250, 0.0980];
green = [0.0784, 0.5804, 0.0784];
lightgreen = [0.0431, 0.4000, 0.1373];

barWidth = 0.2;
figure('position',[100,100,600,500])
b = bar(xdata(2:5:end),[bar_mean2;bar_high2],'stacked','BarWidth', barWidth);
b(1).FaceColor = lightblue;
b(2).FaceColor = blue;
hold on
b = bar(xdata(3:5:end),[bar_mean;bar_high],'stacked','BarWidth', barWidth);
b(1).FaceColor = lightred;
b(2).FaceColor = red;
hold on
b = bar(xdata(4:5:end),[bar_mean_nn;bar_high_nn],'stacked','BarWidth', barWidth); % NN data
b(1).FaceColor = lightgreen;
b(2).FaceColor = green;
hold on
bar(xdata(1:5:end),[bar_mean-bar_mean;bar_high-bar_high],'stacked') 
hold on
bar(xdata(5:5:end),[bar_mean-bar_mean;bar_high-bar_high],'stacked') 
hold on
% bar(xdata(9:5:end),[bar_mean-bar_mean;bar_high-bar_high],'stacked') % NN data
% hold on
grid on
xticks(xdata(3:5:end))
myxtick = [sprintf('% 2d',PN_size_range(1))];
for ii=1:length(PN_size_range)-1
    if PN_size_range(ii+1)<10
        myxtick = [myxtick;sprintf('% 2d',PN_size_range(ii+1))];
    else
        myxtick = [myxtick;sprintf('%2d',PN_size_range(ii+1))];
    end
end
set(gca,'xticklabel',myxtick)
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
ylim([0,30]) 
yticks([2,4,6,8,10,20,30])
xlim([1,max(xdata)-1])
legend('Vanilla RSS-MP [14], (50% percentile)',...
       'Vanilla RSS-MP [14], (90% percentile)',...
       'RSS-MP w/ Dict. Est., (50% percentile)',...
       'RSS-MP w/ Dict. Est., (90% percentile)',...
       'NN Data Driven, (50% percentile)',...
       'NN Data Driven, (90% percentile)')

%% Statistic of training size
figure
plot(dictionary_train_cnt)
xlabel('DFT Beam Indices')
ylabel('Num. as Best Beam')
grid on
