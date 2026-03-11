%% Waterfall Video Reconstruction for Each 5-Frame Set (2D DCT Frame-wise)
clc; clear; close all;

%% --- Load extracted waterfall frames ---
load('video_frames_256x256.mat');   % contains variable 'videoSets'
% videoSets: [256 x 256 x 3 x 5 x 10]
disp('Loaded video frames.');

%% --- Convert to grayscale hyperspectral cube ---
videoSetsGray = zeros(256,256,5,10);
for s = 1:10
    for f = 1:5
        videoSetsGray(:,:,f,s) = rgb2gray(videoSets(:,:,:,f,s));
    end
end

%% --- CASSI Parameters ---
Nx = 256; Ny = 256; L = 5;
p = Nx * Ny * L;
n = 2 * Nx * Ny;
lambda = 1005;
rng(1);

%% --- Function Handles for 2D DCT Operations (Frame-wise) ---
% Forward 2D DCT for each frame: x_frame -> theta_frame
dct2_forward = @(x_frame) dct2_mex(x_frame, Nx, Ny);

% Inverse 2D DCT for each frame: theta_frame -> x_frame  
dct2_inverse = @(theta_frame) idct2_mex(theta_frame, Nx, Ny);

% Combined operators for entire cube
dct_cube_forward = @(x) dct_cube_forward_mex(x, dct2_forward, Nx, Ny, L);
dct_cube_inverse = @(theta) dct_cube_inverse_mex(theta, dct2_inverse, Nx, Ny, L);

%% --- CASSI Measurement Operator ---
% Create random mask for each spectral band
masks = zeros(Nx, Ny, L);
for l = 1:L
    masks(:,:,l) = randi([0,1], Nx, Ny) * 2 - 1;
end

% Function handle for A operator: theta -> y
A_operator = @(theta) cassi_forward_operator_2ddct(theta, masks, dct2_inverse, Nx, Ny, L);

% Function handle for A' operator: y -> theta
A_adjoint = @(y) cassi_adjoint_operator_2ddct(y, masks, dct2_forward, Nx, Ny, L);

%% --- Optimized FISTA with Function Handles ---
fista_function_handle = @(y, lambda, max_iter) ...
    fista_lasso_function_handle_2ddct(y, A_operator, A_adjoint, lambda, max_iter, n, p);

%% --- Loop over all 10 sets ---
for setIdx = 1:10
    fprintf('\n=== Processing Set %d ===\n', setIdx);
    
    % True hyperspectral cube
    X_true = videoSetsGray(:,:,:,setIdx);
    x_true = X_true(:);
    
    % Sparse representation in 2D DCT domain (frame-wise)
    theta_true = dct_cube_forward(x_true);

    %% --- Measurements ---
    y = A_operator(theta_true);
    sigma = 0.01;
    y = y + sigma * randn(n,1);
    
    %% --- FISTA Reconstruction with Function Handles ---
    max_iter = 200;
    theta_hat = fista_function_handle(y, lambda, max_iter);

    %% --- Efficient Debiasing ---
    theta_d = efficient_debiasing_2ddct(y, theta_hat, A_operator, A_adjoint, n, p, sigma);
    
    % Statistical testing (simplified for efficiency)
    alpha = 0.05/p;
    z = norminv(1 - alpha/2);
    B  = @(theta) A_operator(theta);    % = AΨ
    Bt = @(y)     A_adjoint(y);      % = ΨᵀAᵀ

    % The composite operator Q = Psi * Phi * Phi' * Psi'
    Q = @(x) Bt(B(x));   % = (AΨ)'(AΨ x)
    % Variance estimate for θ
    for k=1:1:10
        v = randn(p,1);
        Var_theta = Var_theta+sigma^2 * (1-mu)^2 *Q(v);
    end
    Var_theta = Var_theta/10;

    Var_f = diag(Var_theta);
    Stat = abs(theta_d) * sqrt(n) ./ sigma * sqrt(Var_f);  % Simplified variance estimate
    D = (Stat > z);
    [Sensitivity, Specificity] = compute_sens_spec_new(theta_true, D);

    %% --- Reconstructed cubes ---
    x_hat = dct_cube_inverse(theta_hat);
    x_hat_deb = dct_cube_inverse(theta_d);
    X_hat = reshape(x_hat, Nx, Ny, L);
    X_hat_deb = reshape(x_hat_deb, Nx, Ny, L);

    %% --- Save reconstructed cubes ---
    save(sprintf('reconstruction_set_%02d.mat', setIdx), ...
         'X_true', 'X_hat', 'X_hat_deb', 'theta_hat', 'theta_d', ...
         'Sensitivity', 'Specificity', '-v7.3');

    %% --- Compute and display metrics ---
    snr_fista = compute_snr(X_true, X_hat);
    snr_deb = compute_snr(X_true, X_hat_deb);
    fprintf('  FISTA SNR: %.2f dB, Debiased SNR: %.2f dB\n', snr_fista, snr_deb);
    fprintf('  Sensitivity: %.4f, Specificity: %.4f\n', Sensitivity, Specificity);

    %% --- Visualization ---
    if setIdx <= 3  % Only visualize first 3 sets to save time
        figure('Name', sprintf('Set %d Reconstruction', setIdx), 'NumberTitle','off', ...
               'Position', [100 100 1400 600]);
        
        % Display original frames
        for l = 1:L
            subplot(3, L, l);
            imshow(X_true(:,:,l), []);
            title(sprintf('True Frame %d', l));
        end
        
        % Display FISTA reconstructed frames
        for l = 1:L
            subplot(3, L, L + l);
            imshow(X_hat(:,:,l), []);
            title(sprintf('FISTA Frame %d', l));
        end
        
        % Display debiased reconstructed frames
        for l = 1:L
            subplot(3, L, 2*L + l);
            imshow(X_hat_deb(:,:,l), []);
            title(sprintf('Debiased Frame %d', l));
        end
        
        sgtitle(sprintf('Set %d: True vs FISTA vs Debiased (SNR: %.1f/%.1f dB)', ...
                        setIdx, snr_fista, snr_deb));
    end
end

disp('✅ Completed all 10 sets of 5-frame experiments.');

%% === SUPPORTING FUNCTIONS (2D DCT Frame-wise) ===

function y = cassi_forward_operator_2ddct(theta, masks, dct2_inverse, Nx, Ny, L)
    % Efficient CASSI forward operator using 2D DCT frame-wise
    n = Nx * Ny;
    
    % Convert theta to spatial domain frame by frame
    x = zeros(n * L, 1);
    for l = 1:L
        frame_idx = ((l-1)*n+1):(l*n);
        theta_frame = theta(frame_idx);
        x_frame = dct2_inverse(theta_frame);
        x(frame_idx) = x_frame;
    end
    
    % Apply CASSI forward model
    y = zeros(n, 1);
    for l = 1:L
        frame_idx = ((l-1)*n+1):(l*n);
        x_frame = reshape(x(frame_idx), Nx, Ny);
        masked_frame = x_frame .* masks(:,:,l);
        
        % Apply dispersion (shift)
        shift_amount = l - 1;
        shifted_frame = circshift(masked_frame, [0, shift_amount]);
        y = y + shifted_frame(:);
    end
end

function theta = cassi_adjoint_operator_2ddct(y, masks, dct2_forward, Nx, Ny, L)
    % Efficient CASSI adjoint operator using 2D DCT frame-wise
    n = Nx * Ny;
    y_mat = reshape(y, Nx, Ny);
    
    x_adj = zeros(n * L, 1);
    for l = 1:L
        % Reverse dispersion (shift)
        shift_amount = l - 1;
        shifted_y = circshift(y_mat, [0, -shift_amount]);
        
        % Apply mask
        masked_frame = shifted_y .* masks(:,:,l);
        
        frame_idx = ((l-1)*n+1):(l*n);
        x_adj(frame_idx) = masked_frame(:);
    end
    
    % Convert to 2D DCT domain frame by frame
    theta = zeros(n * L, 1);
    for l = 1:L
        frame_idx = ((l-1)*n+1):(l*n);
        x_frame = x_adj(frame_idx);
        theta_frame = dct2_forward(x_frame);
        theta(frame_idx) = theta_frame;
    end
end

function theta_hat = fista_lasso_function_handle_2ddct(y, A_func, At_func, lambda, max_iter, n, p)
    % FISTA with function handles for 2D DCT frame-wise
    theta = zeros(p, 1);
    theta_prev = theta;
    t = 1;
    L = 0.1;  % Lipschitz constant
    
    fprintf('  FISTA Progress: ');
    
    for iter = 1:max_iter
        % Gradient step using function handles
        residual = A_func(theta) - y;
        grad = At_func(residual);
        theta_temp = theta - (1/L) * grad;
        
        % Soft thresholding
        theta_next = sign(theta_temp) .* max(abs(theta_temp) - lambda/L, 0);
        
        % FISTA acceleration
        t_next = (1 + sqrt(1 + 4*t^2)) / 2;
        theta = theta_next + ((t - 1) / t_next) * (theta_next - theta_prev);
        
        % Update variables
        theta_prev = theta_next;
        t = t_next;
        
        if mod(iter, 50) == 0
            fprintf('%d ', iter);
        end
    end
    fprintf('\n');
    theta_hat = theta_prev;
end

function theta_d = efficient_debiasing_2ddct(y, theta_hat, A_func, At_func, n, p, sigma)
    % Simplified debiasing for 2D DCT frame-wise
    mu = 0.0394;
    
    % Approximate debiasing step
    residual = y - A_func(theta_hat);
    correction = (1/n) * (1-mu) * dct_cube_forward(At_func(residual));
    theta_d = theta_hat + correction;
    
    % Soft threshold to remove very small values
    threshold = sigma * sqrt(2*log(p)) / sqrt(n);
    theta_d(abs(theta_d) < threshold) = 0;
end

function theta = dct2_mex(x_frame, Nx, Ny)
    % 2D DCT for a single frame
    X_frame = reshape(x_frame, Nx, Ny);
    Theta_frame = dct2(X_frame);
    theta = Theta_frame(:);
end

function x = idct2_mex(theta_frame, Nx, Ny)
    % 2D inverse DCT for a single frame
    Theta_frame = reshape(theta_frame, Nx, Ny);
    X_frame = idct2(Theta_frame);
    x = X_frame(:);
end

function theta_cube = dct_cube_forward_mex(x, dct2_forward, Nx, Ny, L)
    % Apply 2D DCT to each frame in the cube
    n = Nx * Ny;
    theta_cube = zeros(n * L, 1);
    
    for l = 1:L
        frame_idx = ((l-1)*n+1):(l*n);
        x_frame = x(frame_idx);
        theta_frame = dct2_forward(x_frame);
        theta_cube(frame_idx) = theta_frame;
    end
end

function x_cube = dct_cube_inverse_mex(theta, dct2_inverse, Nx, Ny, L)
    % Apply 2D inverse DCT to each frame in the cube
    n = Nx * Ny;
    x_cube = zeros(n * L, 1);
    
    for l = 1:L
        frame_idx = ((l-1)*n+1):(l*n);
        theta_frame = theta(frame_idx);
        x_frame = dct2_inverse(theta_frame);
        x_cube(frame_idx) = x_frame;
    end
end

function snr_val = compute_snr(X_true, X_recon)
    % Compute SNR in dB between two hyperspectral cubes
    signal_power = norm(X_true(:))^2;
    noise_power = norm(X_true(:) - X_recon(:))^2;
    if noise_power == 0
        snr_val = inf;
    else
        snr_val = 10 * log10(signal_power / noise_power);
    end
end

function [Sensitivity, Specificity] = compute_sens_spec_new(theta_true, D)
    % Compute sensitivity and specificity for support detection
    true_support = abs(theta_true) > 0;
    detected_support = D > 0;
    
    TP = sum(true_support & detected_support);
    FP = sum(~true_support & detected_support);
    TN = sum(~true_support & ~detected_support);
    FN = sum(true_support & ~detected_support);
    
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
    
    % Handle edge cases
    if (TP + FN) == 0, Sensitivity = 1; end
    if (TN + FP) == 0, Specificity = 1; end
end