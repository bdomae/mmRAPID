clear;clc;
N = 6;
d = 0.57;

PN_beam = (randi(2,1,N)*2-3) + 1j*(randi(2,1,N)*2-3);

total_angle_steps = 999;
angle = linspace(-60,60,total_angle_steps);
angle_stepsize = angle(2) - angle(1);

error_range_len = 6;
MCtimes  = 1e3;
for MCidx = 1:1
    
    % True AoA
    AoA_idx_grid = randi(600) + 200;
    AoA_deg = angle(AoA_idx_grid);
    
    % Steering vector with infinite resolution
    DFT_phase = (0:N-1)*2*pi*d*sind(AoA_deg);
    DFT_beam = exp(-1j * DFT_phase);   
    
    % Region of mainlobe and silobes
    center_idx = AoA_idx_grid;
    near_center_idx = center_idx-66:center_idx+66;
    sidelob_idx = [1:center_idx-67, (center_idx+67):total_angle_steps];
    
    % Array response and pattern without quan
    array_res = exp(1j*2*pi*d*(0:N-1).'*sind(angle))/sqrt(N);
    DFT_pattern = abs(DFT_beam*array_res);
    
    for error_idx = 1:error_range_len
        
        
        
        % Array response and pattern without quan
        quan_scaler = 2*pi/(2^error_idx);
        DFT_phase_quan = round(DFT_phase/quan_scaler) * quan_scaler;
        DFT_beam_quan = exp(-1j * DFT_phase_quan);
        DFT_error_pattern = abs(DFT_beam_quan * array_res);
        
        % Pattner error evaluation
        DFT_gain_loss(MCidx,error_idx) = norm(DFT_pattern(center_idx) - DFT_error_pattern(center_idx))^2/DFT_pattern(center_idx)^2;
        DFT_mainlobe_loss(MCidx,error_idx) = norm(DFT_error_pattern(near_center_idx)...
                                        - DFT_pattern(near_center_idx))^2/norm(DFT_pattern(near_center_idx))^2;
        DFT_sidelobe_loss(MCidx,error_idx) = norm(DFT_error_pattern(sidelob_idx)...
                                        - DFT_pattern(sidelob_idx))^2/norm(DFT_pattern(sidelob_idx))^2;
                                    
        if MCidx == 1 && error_idx == 3
        figure;
        plot(angle(near_center_idx),20*log10(DFT_pattern(near_center_idx)),'k','linewidth',4);
        hold on;
        plot(angle,20*log10(DFT_pattern),'linewidth',2);
        grid on;
        plot(angle,20*log10(DFT_error_pattern),'linewidth',2);
        ylim([-30,30])
        xlabel('Angle [deg]')
        ylabel('Pattern [dB]')
        set(gca,'FontSize',14)
        legend('DFT Beam (mainlobe)','DFT Beam', 'DFT Beam Quan.')
        dim = [.15 .7 .35 .15];
        textcontent = ['Main Cntr Dist ', num2str(10*log10(DFT_gain_loss(MCidx,error_idx)),'%.2f'),...
                       ' [dB], Mainlobe ', num2str(10*log10(DFT_mainlobe_loss(MCidx,error_idx)),'%.2f'),...
                       ' [dB], Sidelobe', num2str(10*log10(DFT_sidelobe_loss(MCidx,error_idx)),'%.2f'),...
                       ' [dB]']
        annotation('textbox',dim,'String',textcontent)
        
        end
    end
    
    
end

% for error_idx = 1:error_range_len
%     error_sigma = error_range(error_idx);
%     DFT_main_loss_bound(error_idx) = (1-cos(error_sigma)^2);
% end

%%
red =      	[0.8500, 0.3250, 0.0980];
figure

plot(1:error_range_len, 10*log10(mean(DFT_gain_loss,1)),'-x','linewidth',2)
hold on
plot(1:error_range_len, 10*log10(mean(DFT_mainlobe_loss,1)),'-o','linewidth',2)
hold on
plot(1:error_range_len, 10*log10(mean(DFT_sidelobe_loss,1)),'-s','linewidth',2)
hold on
grid on
ylim([-60,0])
set(gca,'FontSize',14)
xlabel('Phase Shifter Bits')
ylabel('Distortion Metric [dB]')
lgd = legend('DFT Mainlobe Center (E_0)',...
             'DFT Mainlobe (E_1)',...
             'DFT Sidelobes (E_1)');

lgd.FontSize = 10;