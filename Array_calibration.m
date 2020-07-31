% Alignment performance evaluation

clear;clc;

% parameters
rand('seed',1)

% SNR_dB_pool = -15:5:15;
% SNR_pool = 10.^(SNR_dB_pool./10);
% noise_pow_dB_pool = 0;

MCtimes = 1e1;          % number of monte carlo runs
M = 30;                 % number of actual measurement (noncoherent)
MCS = 10;               % effective number of measurement in complex CS problem
Nr = 36;                % number of antenna snumber
SNR_num = 3;            % number of SNR sweep
SNR_range = linspace(10,20,SNR_num);
HI_flag = 1;            % whether to simulate array phase offset (1 yes, 0 no)
Ref_flag = 1;           % wehther to use Han's approach to refine dictionary in HI
%% Dictionary generation
angle_cand_num = 63;        % candidate angle number in the dictionary
pRX_dict_norm = zeros(angle_cand_num,1);

angle_grid_left = -pi*50/180;
angle_grid_right = pi*50/180;

% Angle grids for the dictionary
angle_grid = linspace(angle_grid_left, angle_grid_right, angle_cand_num);
AOAstep = angle_grid(2) - angle_grid(1);

% Array response dictionary matrix
array_response_dict = exp(1j*(0:Nr-1)'*pi*sin(angle_grid)); 

%% MC simulations

for MCidx = 1:MCtimes
    clc;fprintf('Ite %d out of %d\n',MCidx,MCtimes);
    
    % Regenerate PN codebook to see average effect
    if mod(MCidx,2)==1
        
        % PN sounding codebook
%         W_mat_small = ((randi(2,Nr,MCS)*2-3)+1j*(randi(2,Nr,MCS)*2-3))/sqrt(2*Nr);
%         W_mat = [W_mat_small,W_mat_small,W_mat_small];
        W_mat = ((randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3))/sqrt(2*Nr);

        % Introduces hardware impairment in PN codebook
        if HI_flag
            % add (residual) phase offset in the array
            max_phase_offset = 18/180*pi; % let's say +/- 18 deg
            HI_e = exp(1j * (rand(Nr,1)*2-1) * max_phase_offset);
            W_mat_actual = diag(HI_e) * W_mat;
        else
            % assuming all phase offsets are compensated in calibration
            W_mat_actual = W_mat;
        end
        
        % SVD to decompose APR and ACS
        % Note that you do not know W_vec_rx_actual when HI occurs
        [Ut, Sigmat, Vt] = svd(W_mat.');
        ACS = sqrt(Sigmat(1:MCS,1:MCS)) * Vt(:,1:MCS)';
        
        if HI_flag == 0
            % with perfect calibrated array, APR follows [Yi et al]
            APR = Ut(:,1:MCS) * sqrt(Sigmat(1:MCS,1:MCS));
        else
            if Ref_flag == 0
                % without effort to fix HI, APR follows [Yi et al]
                APR = Ut(:,1:MCS) * sqrt(Sigmat(1:MCS,1:MCS));
            else
                APR = Ut(:,1:MCS) * sqrt(Sigmat(1:MCS,1:MCS));
                
                % Han's heuristic approach to use training data to refine APR
                
                % It requires to collect non-coherent measurement data at
                % all angles probed by W_mat_actual
                y_data = abs(array_response_dict.' * W_mat_actual).^2;
                
                % testtest
                uses_cal = 10;
                BigG = [];
                for rr=1:uses_cal
                    BigG = [BigG; W_mat.'*diag(array_response_dict(:,rr))];
                end
                y_target = reshape(y_data(1:uses_cal,:).',uses_cal*M,1);
                
                cvx_begin sdp quiet
                    variable EE(Nr,Nr) hermitian semidefinite 
                    minimize 0.001*trace(EE)+ norm(diag(BigG * EE * BigG') - y_target, 2)
                cvx_end
                [Umat, Sigma, Vmat] = svd(EE);
                HI_e_hat = sqrt(Sigma(1,1))*Vmat(:,1);
                figure;
                plot(abs(diag(BigG*HI_e_hat*HI_e_hat'*BigG') - y_target))
                
                ACS = sqrt(Sigmat(1:MCS,1:MCS)) * Vt(:,1:MCS)' * diag(HI_e);
                
                
                
%                 temp = Ut(:,1:MCS) * Sigmat(1:MCS,1:MCS)*...
%                        Vt(:,1:MCS)' * diag(HI_e_hat);
%                 norm(temp - W_mat_actual.','fro')/norm(W_mat_actual.','fro')
                
                [U_genie, Sigma_genie, V_genie] = svd(W_mat_actual.');
                ACS_genie = sqrt(Sigma_genie(1:MCS,1:MCS)) * V_genie(:,1:MCS)';
                ACS = ACS_genie;
%                 temp2 = U_genie(:,1:MCS) * Sigma_genie(1:MCS,1:MCS)*...
%                         V_genie(:,1:MCS)';
%                 norm(temp2 - W_mat.','fro')/norm(W_mat.','fro')
                
%                 % Main refinement steo
%                 % Required knowledge:
%                 % 1) ACS (based on knowlegde of W_mat instead of W_mat_actual )
%                 % 2) array_response_dict (array information)  
%                 % 3) y_train (this can be collected from training data)
%                 % Note that it does not require explict info of W_mat_actual
%                 APR = SPR_dictionary_refinement(ACS,...
%                                                 array_response_dict,...
%                                                 y_data);
            end
        end
        
        % Dictionary in RSS matching pursuit (Genie knowledge of HI)
        dict_pRx = abs(W_mat_actual.'*array_response_dict);
        for cc=1:angle_cand_num
            pRX_dict_norm(cc) = norm(dict_pRx(:,cc),2);
        end
        
        % Dictionary in complex CS problem when PhaseLift is used
        complex_CS_dict = ACS * array_response_dict;
        for cc=1:angle_cand_num
            complex_CS_dict_norm(cc) = norm(complex_CS_dict(:,cc),2)^2;
        end
    end
    
    % Common AWGN for all SNR sweep (normalized)
    noise_normal = (randn(M,1)+1j*randn(M,1))/sqrt(2);
    
    for SNRindex = 1:SNR_num
        
        % Noise power in this SNR sweep
        noise_pow = 10^(-SNR_range(SNRindex)/10);
        
        % use a simple channel (LOS at boresight) to verify
        channel = array_response_dict(:,(angle_cand_num+1)/2);
        
        % noncoherent measurement
        y = W_mat_actual.' * channel + noise_normal * sqrt(noise_pow);
        p_vec = abs(y).^2;
        
        % Main solver for SPR beam training
        cvx_begin sdp quiet
            variable Z(MCS,MCS) hermitian semidefinite 
            minimize 0.5*trace(Z) + norm(diag(APR * Z * APR') - p_vec, 2)
        cvx_end
        
        % Takes SVD of Z
        [Umat, Sigma, Vmat] = svd(Z);
        
        % PhaseLift Solution
        cand_score_PLCS = (complex_CS_dict'*Vmat(:,1))./complex_CS_dict_norm';
        [~,bestindex_PLCS(MCidx)] = max(abs(cand_score_PLCS));
        bestAOA_PLCS(MCidx,SNRindex) = angle_grid(bestindex_PLCS(MCidx));
        AOA_error_PLCS(MCidx,SNRindex) = abs(bestAOA_PLCS(MCidx,SNRindex) - 0);
        
        % pRx Matching Pursuit Solution
        cand_score_MP = dict_pRx'*sqrt(p_vec)./pRX_dict_norm;
        [~,bestindex_MP(MCidx)] = max(abs(cand_score_MP));
        bestAOA_MP(MCidx,SNRindex) = angle_grid(bestindex_MP(MCidx));
        AOA_error_MP(MCidx,SNRindex) = abs(bestAOA_MP(MCidx,SNRindex) - 0);
        
    end
end
%% Evaluation of performance

% Evaluate alignment rate
AOAalign_PLCS_mean = sum((AOA_error_PLCS/pi*180)<(105/Nr),1)/MCtimes;
AOAalign_MP_mean = sum((AOA_error_MP/pi*180)<(105/Nr),1)/MCtimes;

figure
plot(SNR_range,AOAalign_PLCS_mean,'-','linewidth',2);
hold on
plot(SNR_range,AOAalign_MP_mean,'-','linewidth',2);
grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
legend('SPR','pRx MP')
set(gca,'FontSize',14)
%% Proposed sparse phase retrieval dictionary refinement from training data and labels
function APR = SPR_dictionary_refinement(ACS, array_response_dict, y_data)
    
    % Estimate dictionary APR using PhaseLift
    % See group meeting slides (Apr. 29) for more details
    
    M = size(y_data,2);
    MCS = M/2;
    for mm=1:M
        % y^H = PN^H(:,mm) * A_rx^H
        % y = A_rx * ACS^H * APR^H(:,mm)
        % where PN = APR * ACS
        
        y_target = y_data(:,mm);

        AthetACS = array_response_dict.' * ACS.';
        cvx_begin sdp quiet
            variable ZZ(MCS,MCS) hermitian semidefinite 
            minimize 1*trace(ZZ)+ norm(diag(AthetACS * ZZ * AthetACS') - y_target, 2)
        cvx_end
        [Umat,Sigma,Vmat] = svd(ZZ);
        APR_hat(:,mm) = Vmat(:,1);

    end
    APR = APR_hat.';
end
