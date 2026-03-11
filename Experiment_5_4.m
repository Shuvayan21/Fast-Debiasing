% Parameters
p = 500;
n_values = [200, 350, 500];
num_trials = 100;
sigma = 1;
xi = 0.1; % correlation coefficient for equicorrelated

% Define Sigma for 3 cases
% (1) Uncorrelated
Sigma_uncorr = sigma^2 * eye(p);

% (2) Band matrix from Javanmard & Montanari (2014), bandwidth = 3
bandwidth = 5;
% Sigma_band = zeros(p);
Sigma_band = sigma^2 * makeSigma(p);

% (3) Equicorrelated
Sigma_equicorr = sigma^2 * ((1-xi) * eye(p) + xi * ones(p));

% Store results
results = cell(3, length(n_values));

% Loop over covariance structures
for cov_case = 1:3
    switch cov_case
        case 1
            Sigma = Sigma_uncorr;
        case 2
            Sigma = Sigma_band;
        case 3
            Sigma = Sigma_equicorr;
    end
    
    % Cholesky factor for generating rows ~ N(0, Sigma)
    L = chol(Sigma, 'lower');
    
    % Loop over n
    for ni = 1:length(n_values)
        n = n_values(ni);
        rho_vals = zeros(num_trials, 1);
        A=zeros(n,p);
        for trial = 1:num_trials
            % Generate A: each row ~ N(0, Sigma)
            for i=1:1:n
                Z = randn(1, p);        % iid N(0,1)
                A(i,:)=Z * L';             % impose covariance Sigma
            end

            % Column norms
            col_norms = sum(A.^2, 1);
            
            % Gram matrix
            G = A' * A;
            G(1:p+1:end) = 0;       % zero out diagonal
            
            % Compute rho(A)
            denom = repmat(col_norms, p, 1); % denominator for each j
            ratios = abs(G) ./ denom;        % elementwise division
            rho = max(ratios(:));            % global max
            
            % Store rho/(1+rho)
            rho_vals(trial) = rho / (1 + rho);
        end
        
        results{cov_case, ni} = rho_vals;
    end
end

%% Plotting with reference lines
figure;
for cov_case = 1:3
    for ni = 1:3
        subplot(3,3,(cov_case-1)*3 + ni);
        histogram(results{cov_case,ni}, 10, 'Normalization', 'probability');
        xlim([0.1 1.1]);   % x-axis 0 to 1
        ylim([0 0.4]);   % y-axis 0 to 1
        hold on;
        
        % Compute threshold
        n = n_values(ni);
        if cov_case == 1
            x_star = 2*sqrt(2) * sqrt(log(p)/n)/(1-2*sqrt(2) * sqrt(log(p)/n));
        elseif cov_case == 2
            Sigma=Sigma_band;
            max_ratio = zeros(p,1);
            for j = 1:p
            % take absolute correlations with column j
            vals = abs(Sigma(:,j)) / Sigma(j,j);
            % remove the diagonal element (l = j)
            vals(j) = -inf;  
            % take maximum
            max_ratio(j) = max(vals);
            end    
            x_star = (2*sqrt(2)+max(max_ratio)) * sqrt(log(p)/n)/(1-2*sqrt(2) * sqrt(log(p)/n));
        else 
            Sigma=Sigma_equicorr;
            max_ratio = zeros(p,1);
            for j = 1:p
            % take absolute correlations with column j
            vals = abs(Sigma(:,j)) / Sigma(j,j);
            % remove the diagonal element (l = j)
            vals(j) = -inf;  
            % take maximum
            max_ratio(j) = max(vals);
            end    
            x_star = (2*sqrt(2)+max(max_ratio)) * sqrt(log(p)/n)/(1-2*sqrt(2) * sqrt(log(p)/n));
        end
        
        % Draw vertical line
        xline(x_star, 'r-', 'LineWidth', 2);
        

        % Title
        switch cov_case
            case 1
                title_str = sprintf('\Sigma_1, n=%d', n);
            case 2
                title_str = sprintf('\Sigma_2, n=%d', n);
            case 3
                title_str = sprintf('\Sigma_3, n=%d', n);
        end
        title(title_str, 'FontSize', 12, 'FontWeight','bold');
         if cov_case==3
            xlabel('\rho/(1+\rho)');
        end
        if ni==1
            ylabel('density');
        end
        set(gca,'FontSize',10,'FontWeight','bold');
        hold off;
    end
end

