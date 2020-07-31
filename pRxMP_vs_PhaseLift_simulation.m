% Alignment performance evaluation

clear;clc;

% parameters
rand('seed',1)

% SNR_dB_pool = -15:5:15;
% SNR_pool = 10.^(SNR_dB_pool./10);
% noise_pow_dB_pool = 0;

runtimes = 3e1;
M = 30;                 % number of measurement
MCS = 15;               % effective number of measurement in complex CS problem
Nr = 36;
SNR_num = 7;
SNR_range = linspace(0,30,SNR_num);

%% dictionary generation
cand_num_r = 63;
dict_pRx_norm = zeros(cand_num_r,1);

cand_y = zeros(M,cand_num_r);
cand_angle_r = linspace(-pi*50/180,pi*50/180,cand_num_r);
AOAstep = cand_angle_r(2) - cand_angle_r(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));

%% MC simulations

for runindex=1:runtimes
    clc;fprintf('Ite %d out of %d\n',runindex,runtimes);
    if mod(runindex,10)==1
        % dictionary generation
        steer_vec_rx = ((randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3))/sqrt(2*Nr);
        steer_vec_rx_actual = diag(exp(1j*2*pi*randn(Nr,1)*0))*steer_vec_rx;
        
        [Ut, Sigmat, Vt] = svd(steer_vec_rx_actual.');
%         ACS = (randn(MCS,Nr) + 1j*randn(MCS,Nr))*sqrt(1/2/Nr);   
        ACS = sqrt(Sigmat(1:MCS,1:MCS))*Vt(:,1:MCS)';
        APR = Ut(:,1:MCS) * sqrt(Sigmat(1:MCS,1:MCS));

%         APR = (randn(M,MCS) + 1j*randn(M,MCS))*sqrt(1/2/MCS);
        
%         Atot_decomp = APR*ACS;
        
        new_dict = cand_ARV_r; 
        
        dict_pRx = abs(steer_vec_rx.'*cand_ARV_r);

        Measure_mat_new = ACS*new_dict;
        for cc=1:cand_num_r
            Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
            dict_pRx_norm(cc) = norm(dict_pRx(:,cc),2);
        end
        
%         % Estimate dictionary APR using PhaseLift
%         A_train = new_dict.';
%         for mm=1:M
%             % y^H = PN^H(:,mm) * A_rx^H
%             % y = A_rx * ACS^H * APR^H(:,mm)
%             % where PN = APR * ACS
%             y_train = abs(A_train * steer_vec_rx(:,mm)).^2;
%         
%             AthetACS = A_train * ACS.';
%             cvx_begin sdp quiet
%                 variable ZZ(MCS,MCS) hermitian semidefinite 
%                 minimize 1*trace(ZZ)+ norm(diag(AthetACS * ZZ * AthetACS') - y_train, 2)
%             cvx_end
%             [Umat,Sigma,Vmat] = svd(ZZ);
%             APR_hat(:,mm) = Vmat(:,1);
%             
%         end
%         W_hat = APR_hat.' * ACS;
        
        
%         norm(abs(A_train * W_hat.') - abs(A_train * steer_vec_rx))^2/...
%         norm(abs(A_train * steer_vec_rx))^2
        
%         figure
%         plot(abs(A_train * W_hat(:,1)));
%         hold on
%         plot(abs(A_train * steer_vec_rx(:,1)))
%         hold on
%         grid on
        
    end
    
    noise_normal = (randn(M,1)+1j*randn(M,1))/sqrt(2);
    for SNRindex = 1:SNR_num
        
        noise_pow = 10^(-SNR_range(SNRindex)/10);
%         APR = APR_hat.';
%         [Ut, Sigmat, Vt] = svd(steer_vec_rx.');
%         ACS = sqrt(Sigmat(1:MCS,1:MCS)) * Vt(:,1:MCS)';
%         y = (APR * ACS * new_dict(:,32)) + noise_normal*sqrt(noise_pow);
        y = (steer_vec_rx.' * new_dict(:,32)) + noise_normal*sqrt(noise_pow);

        d_vec = (abs(y)).^2;

        cvx_begin sdp quiet
            variable Z(MCS,MCS) hermitian semidefinite 
            minimize 0.5 * trace(Z) + norm(diag(APR * Z * APR') - d_vec,2)
        cvx_end
        
%         z_new = (ACS*new_dict(:,32));
%         Z_new = z_new*z_new';
%         norm(diag(APR * Z_new * APR') - d_vec);
        
        [Umat,Sigma,Vmat] = svd(Z);
%         alpha = pinv(Vmat(:,1))*z_new;
%         norm(alpha*Vmat(:,1)-z_new)/norm(z_new);
        
        % PhaseLift Solution
        cand_score_PLCS = (Measure_mat_new'*Vmat(:,1))./Measure_mat_new_norm';
        [~,bestindex_PLCS(runindex)] = max(abs(cand_score_PLCS));
        bestAOA_PLCS(runindex,SNRindex) = (bestindex_PLCS(runindex)-1)*AOAstep-50*pi/180;
        AOA_error_PLCS(runindex,SNRindex) = abs(bestAOA_PLCS(runindex,SNRindex) - 0);
        
        % pRx Matching Pursuit Solution
        cand_score_MP = dict_pRx'*sqrt(d_vec)./dict_pRx_norm;
        [~,bestindex_MP(runindex)] = max(abs(cand_score_MP));
        bestAOA_MP(runindex,SNRindex) = (bestindex_MP(runindex)-1)*AOAstep-50*pi/180;
        AOA_error_MP(runindex,SNRindex) = abs(bestAOA_MP(runindex,SNRindex) - 0);
        
    end
end
%%
AOAalign_PLCS_mean = sum((AOA_error_PLCS/pi*180)<(105/Nr),1)/runtimes;
AOAalign_MP_mean = sum((AOA_error_MP/pi*180)<(105/Nr),1)/runtimes;

figure
plot(SNR_range,AOAalign_PLCS_mean,'-','linewidth',2);
hold on
plot(SNR_range,AOAalign_MP_mean,'-','linewidth',2);
grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
legend('SPR','pRx MP')
set(gca,'FontSize',14)
%%
% figure(99)
% subplot(211)
% [b,a] = ecdf(AOA_error_nocomp(:,5)/pi*180);
% plot(a,b);hold on
% [b,a] = ecdf(AOD_error_mag(:,5)/pi*180);
% plot(a,b);hold on
% grid on
% xlabel('Estimation Error [deg]')
% ylabel('CDF')
% legend('AoD, complex alg','AoD, mag alg')
% xlim([0,105*8/Nt])
% 
% subplot(212)
% [b,a] = ecdf(AOA_error_nocomp(:,5)/pi*180);
% plot(a,b);hold on
% grid on
% xlabel('Estimation Error [deg]')
% ylabel('CDF')
% xlim([0,105*8/Nr])
% legend('AoA, complex alg')
