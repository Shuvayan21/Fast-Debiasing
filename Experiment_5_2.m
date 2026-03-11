%% ====================== PARAMETER SETUP ===============================
p = 500;                               % ambient dimension
n_values = [200, 350, 500];            % sample sizes
num_trials = 100;                      % number of Monte Carlo trials
s = 10;                                % sparsity level
f_sig = 0.05;                          % significance fraction

% Generate support and non-zero coefficients of beta
S = randperm(p, s);
beta = zeros(p,1);
beta(S(1:floor(0.4*s))) = 500 + 500*rand(floor(0.4*s), 1);
beta(S(floor(0.4*s)+1:s)) = 500 + 500*rand(s - floor(0.4*s), 1);

xi = 0.1;                              % equicorrelation coefficient
sigma = 1;

%% ======================= COVARIANCE MATRICES ==========================
Sigma_uncorr = sigma^2 * eye(p);                     % Uncorrelated
Sigma_band   = sigma^2 * makeSigma(p);               % Band matrix
Sigma_equicorr = sigma^2 * ((1-xi)*eye(p) + xi*ones(p)); % Equicorrelated

%% ======================= RESULT STORAGE MATRICES ======================
sens_We = zeros(3, length(n_values));
spec_We = zeros(3, length(n_values));
sens_Dsig = zeros(3, length(n_values));
spec_Dsig = zeros(3, length(n_values));

ETV_We = zeros(3, length(n_values));     % Empirical total variance (W_e)
ATV_We = zeros(3, length(n_values));     % Asymptotic total variance (W_e)
ETV_Dsig = zeros(3, length(n_values));   % Empirical total variance (dSigma)
ATV_Dsig = zeros(3, length(n_values));   % Asymptotic total variance (dSigma)

%% ======================= MAIN SIMULATION LOOP =========================
for cov_case = 1:3
    
    % Select covariance matrix
    switch cov_case
        case 1
            Sigma = Sigma_uncorr;
        case 2
            Sigma = Sigma_band;
        case 3
            Sigma = Sigma_equicorr;
    end
    
    % Cholesky factor for sampling A ~ N(0, Sigma)
    L = chol(Sigma, 'lower');
    
    %% Loop over sample sizes
    for ni = 1:length(n_values)
        n = n_values(ni);
        
        % Pre-allocate matrices for storing debiased estimates
        beta_d_We = zeros(p, num_trials);
        beta_d_Sig = zeros(p, num_trials);
        
        % Generate design matrix A (same for all trials)
        A = randn(n,p) * L';
        
        % Compute required matrices
        rho = compute_rho(A);
        W_e = compute_W_e(A);
        Sigma_beta_We = sigma^2/n * (W_e' * W_e);
        
        %% Monte Carlo trials
        for trial = 1:num_trials
            
            % Noise and observation generation
            eta = sigma * randn(n,1);
            y = A * beta + eta;
            
            %% ====================== LASSO ESTIMATION ======================
            cvx_begin quiet
                variable x_l(p)
                minimize( 0.5 * pow_pos(norm(y - A*x_l), 2) + norm(x_l) )
            cvx_end
            
            %% ============= DEBIASING USING W_e ============================
            beta_l = x_l;
            beta_d_We(:,trial) = beta_l + (1/n) * W_e' * (y - A*beta_l);
            
            % Compute t-statistics and classification
            TG = beta_d_We(:,trial) ./ sqrt(diag(Sigma_beta_We));
            J = abs(TG) > 2.33;     % Hypothesis test at 1% level
            
            infec = (abs(beta) > 0);
            
            % Confusion counts
            tp = sum(J==1 & infec==1);
            fp = sum(J==1 & infec==0);
            tn = sum(J==0 & infec==0);
            fn = sum(J==0 & infec==1);
            
            % Store metrics
            sens_We(cov_case,ni) = sens_We(cov_case,ni) + tp/(tp+fn);
            spec_We(cov_case,ni) = spec_We(cov_case,ni) + tn/(tn+fp);
            
            %% ============= DEBIASING USING dSigma^{-1} ====================
            W = compute_d_Sigma(A, Sigma, beta);
            Sigma_beta_W = sigma^2/n * (W' * W);
            
            beta_d_Sig(:,trial) = beta_l + (1/n) * W' * (y - A*beta_l);
            
            % t-statistics
            TG = beta_d_Sig(:,trial) ./ sqrt(diag(Sigma_beta_W));
            J = abs(TG) > 2.33;
            
            % Confusion counts
            tp = sum(J==1 & infec==1);
            fp = sum(J==1 & infec==0);
            tn = sum(J==0 & infec==0);
            fn = sum(J==0 & infec==1);
            
            % Store metrics
            sens_Dsig(cov_case,ni) = sens_Dsig(cov_case,ni) + tp/(tp+fn);
            spec_Dsig(cov_case,ni) = spec_Dsig(cov_case,ni) + tn/(tn+fp);
        end
        
        %% ======================= AVERAGING RESULTS =======================
        sens_We(cov_case,ni)  = sens_We(cov_case,ni) / num_trials;
        spec_We(cov_case,ni)  = spec_We(cov_case,ni) / num_trials;
        sens_Dsig(cov_case,ni)= sens_Dsig(cov_case,ni) / num_trials;
        spec_Dsig(cov_case,ni)= spec_Dsig(cov_case,ni) / num_trials;
        
        % Variance computations
        ETV_We(cov_case,ni)  = mean(var(beta_d_We, 0, 2));
        ETV_Dsig(cov_case,ni)= mean(var(beta_d_Sig, 0, 2));
        ATV_We(cov_case,ni)  = mean(diag(Sigma_beta_We));
        ATV_Dsig(cov_case,ni)= mean(diag(Sigma_beta_W));
        
    end
end

%% ======================= GENERATE RESULT TABLES =========================

% Column names
VarNames = {'n', 'sens_We', 'spec_We', 'sens_Dsig', 'spec_Dsig', 'ETV_Ratio', 'ATV_Ratio'};

% Table for each covariance case
for cov_case = 1:3
    T = table( ...
        n_values', ...
        sens_We(cov_case,:)', ...
        spec_We(cov_case,:)', ...
        sens_Dsig(cov_case,:)', ...
        spec_Dsig(cov_case,:)', ...
        (ETV_We(cov_case,:) ./ ETV_Dsig(cov_case,:))', ...
        (ATV_We(cov_case,:) ./ ATV_Dsig(cov_case,:))', ...
        'VariableNames', VarNames );
    
    fprintf("\n======================== TABLE %d ========================\n", cov_case);
    disp(T)
end