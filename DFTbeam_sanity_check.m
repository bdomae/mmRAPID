clear;clc;

% Concantenate all csv files into one variable
x = [];
for batch_idx = 1:5
    % Need to manually delete column for 'link_dir', 'tx_beam_label', 'rx_beam_label'
    % 'tx_mask', and 'rx_mask' in this csv files
    
    file_name = ['results_awv0_20-04-15_part',num2str(batch_idx),'_ext.csv'];
    x_temp = csvread(file_name,1,0);
    x = [x;x_temp];
end

% Parameters
target_filter_beam_idx = 2; % change this to change DFT peak detection target
PN_beam_num = 36;            % fixed
DFT_beam_num = 64;           % fixed
axis_offset = 2000;          % workaround for plotting; don't change this
data_len = size(x,1);        % length of data
x_rx_beam = x(1:data_len,2); % after modification, the 4-th column should be rx_beam_idx
x_pRx = x(1:data_len,18);    % after modification, the 18-th column should be pRx
cnt = 1;
x_axis_angle = [linspace(-45,45,DFT_beam_num),...
                linspace(axis_offset+1,axis_offset+PN_beam_num,PN_beam_num)];
figure('Position',[100 100 1200 800])
zz=1;
legendtext={}

% DFT peak detection, filter, and plot
for ii=1:data_len-1
    if x_rx_beam(ii)>x_rx_beam(ii+1)
        data_to_plot = [x_axis_angle(1+x_rx_beam(cnt:ii))',x_pRx(cnt:ii)];
        [~,pRx_max_idx] = max(x_pRx(cnt:ii));
        temp = x_rx_beam(cnt:ii);
        if temp(pRx_max_idx) == target_filter_beam_idx
            [~,zzz] = sort(x_axis_angle(1+x_rx_beam(cnt:ii)));
            subplot(211)
            plot(data_to_plot(zzz,1),data_to_plot(zzz,2),'linewidth',2);
            legendtext{zz} = ['position', num2str((zz-1)*3+1,'%d')];
            hold on
            subplot(212)
            plot(data_to_plot(zzz,1)-axis_offset,10.^(data_to_plot(zzz,2)/10),'linewidth',2);
            legendtext{zz} = ['position', num2str((zz-1)*3+1,'%d')];
            hold on
            zz = zz+1;
        end
        cnt = ii;
    end
    if ii+1 == data_len
        data_to_plot = [x_axis_angle(1+x_rx_beam(cnt:ii))',x_pRx(cnt:ii)];
        [~,zzz] = sort(x_axis_angle(1+x_rx_beam(cnt:ii)));
%         plot(data_to_plot(zzz,1),data_to_plot(zzz,2),'linewidth',2);
%         legendtext{zz} = ['position', num2str((zz-1)*3+1,'%d')];
%         hold on
        cnt = ii;
    end
end

% Setting for plot
subplot(211)
xlabel('Angle [deg]')
ylabel('pRx [dB]')
grid on
xlim([-45,45])
legend(legendtext,'Location','eastoutside','FontSize',11,'NumColumns',2)

subplot(212)
xlabel('PN Beam Index')
ylabel('pRx [non-log scale]')
grid on
xlim([1,PN_beam_num])
legend(legendtext,'Location','eastoutside','FontSize',11,'NumColumns',2)