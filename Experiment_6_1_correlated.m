clear; clc; close all;

%% PARAMETERS
n = 20000;                 % number of measurements
p_side = 256;              % image size 256x256
N = p_side;
p = p_side^2;

sigma_factor = 0.01;
K = 50;                    % repeated trials
alpha = 0.05;
lambda = 1e-2;             % tune this

%% DCT forward/inverse operators (matrix-free)
dct2_handle  = @(X) dct2(X);
idct2_handle = @(X) idct2(X);

%% Image loop
img_names = {'barbara','cameraman','moon','male'};

for idx = 1:length(img_names)

    fprintf("\n======= Processing %s =======\n", img_names{idx});

    %% Load and vectorize image
    f = im2double(imresize(imread([img_names{idx} '.png']), [N N]));
    if size(f,3)>1, f = rgb2gray(f); end
    f_vec = f(:);

    %% Compute DCT coefficients θ_true
    theta_true_mat = dct2(f);
    theta_true = theta_true_mat(:);

    %% --- SENSING MATRIX WITH CROSSTALK -----------------------
    % A is Rademacher sensing matrix
    A = generate_crosstalk_matrix(n, p, 0.02);


    %% Operator handles for A_tilde Ψ
    A_func     = @(x) A * x;
    At_func    = @(y) A * y;

    A_Psi_func = @(theta) A_func( reshape(idct2_handle( reshape(theta, N, N) ), [], 1) );
    PsiT_AT_func = @(y) reshape( dct2_handle( reshape(At_func(y), N, N) ), [], 1 );

    %% Generate noisy measurements
    y_clean = A_Psi_func(theta_true(:));
    sigma = sigma_factor * mean(abs(y_clean));
    y = y_clean + sigma * randn(n,1);

    %% FISTA options
    opts.pos = false;
    opts.max_iter = 200;
    opts.tol = 1e-6;
    opts.lambda = lambda;
    opts.backtracking = false;

    %% Run Lasso (FISTA)
    fprintf("Running FISTA Lasso...\n");
    theta_hat = fista_lasso_function_handle(y, A_Psi_func, PsiT_AT_func, n, p, opts);

    %% --- FAST DEBIASING --------------------------------------
    fprintf("Running fast debiasing...\n");

    % Compute Wtilde approx = (AΨ)' / (n * column norms of AΨ)^2
    % Column norms estimated implicitly
    rho = compute_rho_function_handle(A_func, At_func, n, p);
    mu = rho/(1+rho);

    % Debiasing operation as function handle
    debias_func = @(theta) theta + (1/n) * (1-mu) * PsiT_AT_func(y - A_Psi_func(theta));
    theta_deb = debias_func(theta_hat);
    % reconstruct
    f_hat_lasso = idct2(reshape(theta_hat, N, N));
    f_hat_deb   = idct2(reshape(theta_deb, N, N));

    %% Coverage estimation (matrix-free)
    fprintf("Computing coverage...\n");

    zval = norminv(1 - alpha/2);
    B  = @(theta) A_Psi_func(theta);    % = AΨ
    Bt = @(y)     PsiT_AT_func(y);      % = ΨᵀAᵀ

    % The composite operator Q = Psi * Phi * Phi' * Psi'
    Q = @(x) Bt(B(x));   % = (AΨ)'(AΨ x)
    % Variance estimate for θ
    for k=1:1:10
        v = randn(p,1);
        Var_theta = Var_theta+sigma^2 * (1-mu)^2 *Q(v);
    end
    Var_theta = Var_theta/10;

    Var_f = diag(Var_theta);

    CI_half = zval * sqrt(abs(Var_f));

    coverage_count = zeros(p,1);

    for k = 1:K
        eta = sigma * randn(n,1);
        yk = y_clean + eta;

        theta_k = theta_hat + Wtilde_apply(yk - A_Psi_func(theta_hat));
        f_k = idct2(reshape(theta_k, N, N));
        f_k = f_k(:);

        inside = abs(f_vec - f_k) <= CI_half;
        coverage_count = coverage_count + inside;
    end

    coverage_map = coverage_count / K;
    coverage_img = reshape(coverage_map, N, N);

    %% Metrics
    RRMSE = norm(theta_deb(:) - theta_true(:)) / norm(theta_true(:));
    SSIMv = ssim(f_hat_deb, f);

    edges = edge(f, 'Canny');
    edge_locs = find(edges(:));
    non_edge_locs = setdiff(1:p, edge_locs);

    cov_edges = mean(coverage_map(edge_locs));
    cov_non_edges = mean(coverage_map(non_edge_locs));

    fprintf("RRMSE: %.4f\n", RRMSE);
    fprintf("SSIM:  %.4f\n", SSIMv);
    fprintf("Coverage (edges): %.3f\n", cov_edges);
    fprintf("Coverage (non-edges): %.3f\n", cov_non_edges);

    %% Plots
    figure;
    subplot(1,4,1); imshow(f, []); title("Original");
    subplot(1,4,2); imshow(f_hat_lasso, []); title("Lasso");
    subplot(1,4,3); imshow(f_hat_deb, []); title("Debiased");
    subplot(1,4,4); imagesc(coverage_img); axis off; colorbar; title("Coverage");

end 

function theta = solve_lasso(Phi,y,lambda)
    [B,FitInfo] = lasso(Phi,y,'Lambda',lambda,'Standardize',false);
    theta = B;
end

function D = make_dct_matrix(n)
    D = kron(dctmtx(n), dctmtx(n));
end

function Var_f = compute_pixel_variances(Var_theta, N)
    Psi = make_dct_matrix(N);
    Var_f = sigma^2 * sum((Psi.^2) .* Var_theta', 2);
end

function coverage_map = compute_coverage(Phi, Psi, PsiT, f_true, theta_hat, theta_true, sigma, Wtilde, alpha, K, p_side);
    p = length(f_true);
    zval = norminv(1 - alpha/2);
    Var_theta = sigma^2 * diag(Wtilde'*Wtilde);
    Var_f = compute_pixel_variances(Var_theta, sqrt(p));
    CI_half = zval * sqrt(Var_f);

    %f_hat = reshape(PsiT(reshape(theta_deb,N,N)), p, 1);

    coverage_count = zeros(p,1);

    for k = 1:K
        eta = sigma * randn(size(Phi,1),1);
        yk = Phi*(Psi(f_true)) + eta;
        % No need to re-run Lasso; use linearized debiasing:
        theta_k = theta_hat + (1/numel(eta))*Wtilde'*(yk - Phi*(Psi(theta_deb)));
        f_hat = reshape(PsiT(reshape(theta_k,N,N)), p,1);

        inside = abs(f_true - f_hat) <= CI_half;
        coverage_count = coverage_count + inside;
    end

    coverage_map = coverage_count / K;
end

function theta_hat = fista_lasso_function_handle(y, A_func, AT_func, n, p, opts)
    % FISTA with function handles to avoid large matrix storage
    max_iter = opts.max_iter;
    lambda = opts.lambda;
    tol = opts.tol;
    
    % Initialize
    theta = zeros(p, 1);
    theta_prev = theta;
    t = 1;
    t_prev = 1;
    
    % Lipschitz constant estimation
    %L = estimate_lipschitz(A_func, AT_func, n, p);
    L=23;
    for iter = 1:max_iter
        % Gradient step
        residual = A_func(theta) - y;
        gradient = AT_func(residual);
        theta_temp = theta - (1/L) * gradient;
        
        % Soft thresholding
        theta_next = sign(theta_temp) .* max(0, abs(theta_temp) - lambda/L);
        
        % FISTA acceleration
        t_next = (1 + sqrt(1 + 4*t^2)) / 2;
        theta = theta_next + ((t - 1) / t_next) * (theta_next - theta_prev);
        
        % Check convergence
        if norm(theta - theta_prev) / norm(theta_prev) < tol
            break;
        end
        
        % Update variables
        theta_prev = theta_next;
        t_prev = t;
        t = t_next;
        
        if mod(iter, 10) == 0
            fprintf('Iteration %d/%d\n', iter, max_iter);
        end
    end
    
    theta_hat = theta_next;
end

function L = estimate_lipschitz(A_func, AT_func, n, p)
    % Power iteration to estimate Lipschitz constant
    x = randn(p, 1);
    x = x / norm(x);
    
    for i = 1:5  % Reduced iterations for speed
        Ax = A_func(x);
        x = AT_func(Ax);
        L_est = norm(x);
        x = x / L_est;
    end
    
    L = L_est * 1.1; % Safety factor
end

function rho = compute_rho_function_handle(A_func, AT_func, n, p)
    % Estimate rho parameter without storing large matrices
    % Using power iteration to find maximum eigenvalue of A'A
    x = randn(p, 1);
    x = x / norm(x);
    
    for i = 1:5  % Reduced iterations for speed
        Ax = A_func(x);
        x = AT_func(Ax);
        sigma_max = norm(x);
        x = x / sigma_max;
    end
    
    rho = sigma_max^2 / n;
end