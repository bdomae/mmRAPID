clear;clc;
N = 36;
d = 0.57;
DFT_phase = (0:N-1)*2*pi*d*sind(0);
% DFT_phase_quan = round(DFT_phase/(pi/4))*(pi/4);
DFT_beam = exp(-1j * DFT_phase);
% DFT_beam_quan = exp(-1j * DFT_phase_quan);
PN_beam = (randi(2,1,N)*2-3) + 1j*(randi(2,1,N)*2-3);
center_idx = 180;
near_center_idx = center_idx-30:center_idx+30;
sidelob_idx = [1:center_idx-3,center_idx+3:center_idx*2-1];
angle = linspace(-45,45,center_idx*2-1);
error_range_len = 20;
error_range = linspace(0.1,20,error_range_len)/180*pi;
MCtimes  = 1e3;
for MCidx = 1:MCtimes
    for error_idx = 1:error_range_len
        error_sigma = error_range(error_idx);
        mismatch(error_idx,:) = exp(1j*(rand(1,N)*2-1)*error_sigma);
%         mismatch(error_idx,:) = ones(1,N);
%         mismatch(error_idx,8) = 2;
%         mismatch(error_idx,20) = 2;
%         mismatch(error_idx,32) = 2;
    end
    array_res = exp(1j*2*pi*d*(0:N-1).'*sind(angle))/sqrt(N);
    
    DFT_error_pattern = abs(mismatch * diag(DFT_beam)*array_res);
    PN_error_pattern = abs(mismatch * diag(PN_beam)*array_res);
    
    DFT_pattern = abs(DFT_beam*array_res);
    PN_pattern = abs(PN_beam*array_res);
    
    % Two-element calibration
    for error_idx = 1:error_range_len
        error_sigma = error_range(error_idx);
        cal_loss(error_idx) = 20*log10(2) - 20*log10(abs(1 + exp(1j*error_sigma)));
    end
    
    for error_idx = 1:error_range_len
        DFT_gain_loss(MCidx,error_idx) = norm(DFT_pattern(center_idx) - DFT_error_pattern(error_idx,center_idx))^2/DFT_pattern(center_idx)^2;
        DFT_mainlobe_loss(MCidx,error_idx) = norm(DFT_error_pattern(error_idx,near_center_idx)...
                                        - DFT_pattern(near_center_idx))^2/norm(DFT_pattern(near_center_idx))^2;
        DFT_sidelobe_loss(MCidx,error_idx) = norm(DFT_error_pattern(error_idx,sidelob_idx)...
                                        - DFT_pattern(sidelob_idx))^2/norm(DFT_pattern(sidelob_idx))^2;
        PN_loss(MCidx,error_idx) = norm(PN_error_pattern(error_idx,:) - PN_pattern)^2/norm(PN_pattern)^2;
    end
end

% for error_idx = 1:error_range_len
%     error_sigma = error_range(error_idx);
%     DFT_main_loss_bound(error_idx) = (1-cos(error_sigma)^2);
% end

%%
red =      	[0.8500, 0.3250, 0.0980];
figure
% subplot(211)
% plot(error_range/pi*180, -cal_loss,'k--','linewidth',2)
% xlabel('Max Element-wise Phase Offset [deg]')
% ylabel('Gain Loss during Calibration [dB]')
% hold on
% grid on
% subplot(212)

plot(error_range/pi*180, 10*log10(mean(DFT_mainlobe_loss,1)),'-o','linewidth',2)
hold on
% plot(error_range/pi*180, 10*log10(mean(DFT_gain_loss,1)),'-o','linewidth',2,'Color',red)
hold on
% plot(linspace(0.1,20,100), 10*log10((1-cosd(linspace(0.1,20,100))).^2),'--','linewidth',2,'Color',red)
% b.s
hold on
plot(error_range/pi*180, 10*log10(mean(DFT_sidelobe_loss,1)),'-s','linewidth',2)
hold on
plot(error_range/pi*180, 10*log10(mean(PN_loss,1)),'-x','linewidth',2)
grid on
ylim([-60,0])
set(gca,'FontSize',14)
xlabel('Max Phase Offset \psi_m_a_x [deg]')
ylabel('Distortion Metric [dB]')
lgd = legend('DFT Mainlobe (E_0)',...
       'DFT Sidelobes (E_1)',...
       'PN All Lobes (E_2)');
% lgd = legend('DFT Mainlobe (E_0)',...
%        'DFT Mainlobe Center (E_1)',...
%        'DFT Mainlobe Center Bound (E_2)',...
%        'DFT Sidelobes (E_3)',...
%        'PN All Lobes (E_4)');
lgd.FontSize = 10;