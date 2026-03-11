%% CASSI with FISTA Lasso reconstruction - RECONSTRUCT 24 SPECTRAL BANDS
clear; close all; clc;

%% Load the data
fprintf('Loading measurement data...\n');
scene_data = load('Y_subset.mat');
y_measurement = double(scene_data.Y_subset); 

fprintf('Loading mask data...\n');
mask_data = load('Cu_subset.mat');
mask = double(mask_data.Cu_subset);

[H, W, num_frames] = size(y_measurement);
[~, ~, num_spectral_bands, num_temporal_instances] = size(mask);

fprintf('Measurement size: %d x %d x %d (spatial x spatial x temporal frames)\n', H, W, num_frames);
fprintf('Mask size: %d x %d x %d x %d (spatial x spatial x spectral x temporal)\n', H, W, num_spectral_bands, num_temporal_instances);
fprintf('Target reconstruction: %d x %d x %d (spatial x spatial x spectral)\n', H, W, num_spectral_bands);

assert(num_frames == num_temporal_instances, 'Number of measurement frames must match temporal instances in mask');
assert(num_spectral_bands == num_temporal_instances, 'Number of spectral bands must match temporal instances');

%% Use ALL temporal frames for reconstruction
fprintf('Using ALL %d temporal frames to reconstruct %d spectral bands...\n', num_frames, num_spectral_bands);

%% Vectorize ALL measurements for joint reconstruction
fprintf('Vectorizing all %d temporal measurements...\n', num_frames);

% Vectorize all measurements
y_all = [];
for t = 1:num_frames
    y_frame = y_measurement(:,:,t);
    y_all = [y_all; y_frame(:)];
end

M_total = length(y_all);
N = H * W * num_spectral_bands;  % Size of hyperspectral cube to reconstruct

fprintf('Total measurement vector size: %d\n', M_total);
fprintf('Hyperspectral cube vector size: %d\n', N);
fprintf('Compression ratio: %.2f:1\n', N/M_total);

%% Create joint forward and adjoint operators using ALL temporal frames
fprintf('Creating joint forward/adjoint operators using all %d frames...\n', num_frames);

A = @(x) forward_op_joint_spectral(x, mask, H, W, num_spectral_bands, num_frames);
At = @(y) adjoint_op_joint_spectral(y, mask, H, W, num_spectral_bands, num_frames);

%% Verify forward-adjoint consistency
fprintf('Verifying forward-adjoint consistency...\n');
test_x = randn(N, 1);
test_y = A(test_x);
test_At = At(test_y);
consistency = abs(test_x' * test_At - norm(test_y)^2) / (norm(test_x) * norm(test_y));
fprintf('Forward-adjoint consistency error: %.2e\n', consistency);

%% FISTA parameters for joint reconstruction
params = struct();
params.max_iter = 100;
params.lambda = 0.001;  % Regularization parameter
params.tol = 1e-5;
params.verbose = true;

%% Run JOINT FISTA reconstruction using ALL frames
fprintf('Starting JOINT FISTA reconstruction to recover %d spectral bands...\n', num_spectral_bands);
tic;
reconstructed_HSI = fista_lasso_diagonal(y_all, A, At, H, W, num_spectral_bands, params);
total_time = toc;

fprintf('Joint reconstruction completed in %.2f seconds\n', total_time);
fprintf('Successfully reconstructed %d spectral bands using %d coded snapshots\n', num_spectral_bands, num_frames);

%% Display comprehensive results
fprintf('Displaying comprehensive reconstruction results...\n');

% Define wavelengths for 24 bands (adjust range as needed)
wavelengths = linspace(450, 650, num_spectral_bands);

% Display the complete hyperspectral reconstruction
display_spectral_reconstruction_results(reconstructed_HSI, y_measurement, mask, wavelengths);

% Analyze reconstruction quality
analyze_spectral_reconstruction_quality(reconstructed_HSI, y_measurement, mask, wavelengths);

% Save complete results
fprintf('Saving complete spectral reconstruction...\n');
save('spectral_reconstruction_24bands.mat', ...
     'reconstructed_HSI', 'y_measurement', 'mask', 'wavelengths', ...
     'total_time', '-v7.3');

fprintf('Complete spectral reconstruction with %d bands saved!\n', num_spectral_bands);

%% Display results
fprintf('Displaying reconstruction results...\n');

L=24;
    for l=1:L
        wavelength(l) = 450 + (l-1) * (650-450) / (L-1);
    end
    HSI_hat = upscale_hyperspectral_image(reconstructed_HSI);
    plot_colored_cube_single_figure(HSI_hat, 500, wavelength, 1:L)
%% Support functions for spectral reconstruction

function y = forward_op_joint_spectral(x, mask, H, W, num_bands, num_frames)
    % Joint forward operator using ALL temporal frames
    % For each temporal frame t: y_t = sum_{l=1}^{num_bands} diag(mask_{t,l}) * x_l
    % where x_l is the l-th spectral band
    
    % Reshape x to hyperspectral cube [H, W, num_bands]
    X = reshape(x, [H, W, num_bands]);
    y = [];
    
    for t = 1:num_frames
        y_frame = zeros(H, W);
        for l = 1:num_bands
            % Use mask for temporal instance t and spectral band l
            y_frame = y_frame + mask(:, :, l, t) .* X(:, :, l);
        end
        y = [y; y_frame(:)];
    end
end

function x = adjoint_op_joint_spectral(y, mask, H, W, num_bands, num_frames)
    % Joint adjoint operator using ALL temporal frames
    % For each temporal frame t: x_l += diag(mask_{t,l}) * y_t for each band l
    
    x = zeros(H, W, num_bands);
    frame_size = H * W;
    
    for t = 1:num_frames
        % Extract measurement for this temporal frame
        start_idx = (t-1) * frame_size + 1;
        end_idx = t * frame_size;
        y_frame = reshape(y(start_idx:end_idx), [H, W]);
        
        % Accumulate adjoint operation for each spectral band
        for l = 1:num_bands
            x(:, :, l) = x(:, :, l) + mask(:, :, l, t) .* y_frame;
        end
    end
    x = x(:);
end

function display_spectral_reconstruction_results(reconstructed_HSI, y_measurement, mask, wavelengths)
    % Display spectral reconstruction results
    [H, W, num_bands] = size(reconstructed_HSI);
    [~, ~, ~, num_frames] = size(mask);
    
    figure('Position', [50, 50, 1800, 1000]);
    
    % 1. Original encoded measurements (first few frames)
    subplot(3, 4, 1);
    imagesc(y_measurement(:,:,1));
    axis image; colorbar;
    title('Coded Snapshot 1');
    colormap(parula);
    
    subplot(3, 4, 2);
    imagesc(y_measurement(:,:,round(num_frames/2)));
    axis image; colorbar;
    title(sprintf('Coded Snapshot %d', round(num_frames/2)));
    colormap(parula);
    
    subplot(3, 4, 3);
    imagesc(y_measurement(:,:,num_frames));
    axis image; colorbar;
    title(sprintf('Coded Snapshot %d', num_frames));
    colormap(parula);
    
    % 2. Reconstructed spectral bands
    subplot(3, 4, 4);
    rgb_composite = create_spectral_rgb_composite(reconstructed_HSI, wavelengths);
    imshow(rgb_composite);
    title('Reconstructed: RGB Composite');
    
    % 3. Representative spectral bands
    subplot(3, 4, 5);
    % First band (blue)
    imagesc(reconstructed_HSI(:,:,1));
    axis image; colorbar;
    title(sprintf('Band 1: %.0f nm', wavelengths(1)));
    colormap(jet);
    
    subplot(3, 4, 6);
    % Middle band (green)
    mid_band = round(num_bands/2);
    imagesc(reconstructed_HSI(:,:,mid_band));
    axis image; colorbar;
    title(sprintf('Band %d: %.0f nm', mid_band, wavelengths(mid_band)));
    colormap(jet);
    
    subplot(3, 4, 7);
    % Last band (red)
    imagesc(reconstructed_HSI(:,:,num_bands));
    axis image; colorbar;
    title(sprintf('Band %d: %.0f nm', num_bands, wavelengths(num_bands)));
    colormap(jet);
    
    subplot(3, 4, 8);
    % Mean across all spectral bands
    mean_spectral = mean(reconstructed_HSI, 3);
    imagesc(mean_spectral);
    axis image; colorbar;
    title('Mean Across Spectral Bands');
    colormap(parula);
    
    % 4. Spectral profiles at different locations
    subplot(3, 4, 9);
    locations = [round(H/4), round(W/4); 
                 round(H/2), round(W/2);
                 round(3*H/4), round(3*W/4)];
    
    hold on;
    colors = ['r', 'g', 'b'];
    for i = 1:size(locations, 1)
        h = locations(i, 1);
        w = locations(i, 2);
        spectral_profile = squeeze(reconstructed_HSI(h, w, :));
        plot(wavelengths, spectral_profile, colors(i), 'LineWidth', 2);
    end
    xlabel('Wavelength (nm)');
    ylabel('Intensity');
    title('Spectral Profiles');
    legend('Location 1', 'Location 2', 'Location 3');
    grid on;
    
    % 5. Data statistics
    subplot(3, 4, 10);
    histogram(reconstructed_HSI(:), 50);
    xlabel('Intensity');
    ylabel('Frequency');
    title('Intensity Distribution');
    grid on;
    
    % 6. Verify reconstruction with first coded snapshot
    subplot(3, 4, 11);
    simulated_meas = zeros(H, W);
    for l = 1:num_bands
        simulated_meas = simulated_meas + mask(:,:,l,1) .* reconstructed_HSI(:,:,l);
    end
    imagesc(simulated_meas);
    axis image; colorbar;
    title('Simulated Coded Snapshot 1');
    colormap(parula);
    
    % 7. Reconstruction info
    subplot(3, 4, 12);
    text(0.1, 0.9, sprintf('RECONSTRUCTION SUMMARY'), 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.7, sprintf('Spatial size: %d x %d', H, W));
    text(0.1, 0.6, sprintf('Spectral bands: %d', num_bands));
    text(0.1, 0.5, sprintf('Coded snapshots: %d', num_frames));
    text(0.1, 0.4, sprintf('Wavelength range: %.0f-%.0f nm', min(wavelengths), max(wavelengths)));
    text(0.1, 0.3, sprintf('Mean intensity: %.3f', mean(reconstructed_HSI(:))));
    text(0.1, 0.2, sprintf('Dynamic range: [%.3f, %.3f]', min(reconstructed_HSI(:)), max(reconstructed_HSI(:))));
    axis off;
    
    sgtitle(sprintf('Spectral Reconstruction: %d Bands from %d Coded Snapshots', num_bands, num_frames), ...
            'FontSize', 16, 'FontWeight', 'bold');
end

function analyze_spectral_reconstruction_quality(reconstructed_HSI, y_measurement, mask, wavelengths)
    % Analyze reconstruction quality
    [H, W, num_bands] = size(reconstructed_HSI);
    [~, ~, ~, num_frames] = size(mask);
    
    figure('Position', [100, 100, 1400, 600]);
    
    % Calculate residuals for each temporal frame
    residuals = zeros(H, W, num_frames);
    mse_per_frame = zeros(num_frames, 1);
    
    for t = 1:num_frames
        % Simulate measurement from reconstruction
        simulated_meas = zeros(H, W);
        for l = 1:num_bands
            simulated_meas = simulated_meas + mask(:,:,l,t) .* reconstructed_HSI(:,:,l);
        end
        residuals(:,:,t) = y_measurement(:,:,t) - simulated_meas;
        mse_per_frame(t) = mean(mean(residuals(:,:,t).^2));
    end
    
    % Display residuals
    subplot(2, 3, 1);
    imagesc(mean(residuals, 3));
    axis image; colorbar;
    title('Mean Residual Across Frames');
    colormap(jet);
    
    subplot(2, 3, 2);
    plot(1:num_frames, mse_per_frame, 'ro-', 'LineWidth', 2);
    xlabel('Temporal Frame');
    ylabel('MSE');
    title('Reconstruction Error per Frame');
    grid on;
    
    subplot(2, 3, 3);
    histogram(residuals(:), 50);
    xlabel('Residual Value');
    ylabel('Frequency');
    title('Residual Distribution');
    grid on;
    
    % Spectral characteristics
    subplot(2, 3, 4);
    mean_spectrum = squeeze(mean(mean(reconstructed_HSI, 1), 2));
    plot(wavelengths, mean_spectrum, 'b-', 'LineWidth', 2);
    xlabel('Wavelength (nm)');
    ylabel('Mean Intensity');
    title('Mean Spectral Signature');
    grid on;
    
    % Spatial characteristics
    subplot(2, 3, 5);
    spatial_variance = std(reconstructed_HSI, 0, 3);
    imagesc(spatial_variance);
    axis image; colorbar;
    title('Spectral Variance Map');
    colormap(jet);
    
    % Quality metrics
    subplot(2, 3, 6);
    overall_mse = mean(mse_per_frame);
    snr = 20 * log10(max(reconstructed_HSI(:)) / sqrt(overall_mse));
    sparsity_ratio = nnz(reconstructed_HSI) / numel(reconstructed_HSI);
    
    text(0.1, 0.8, sprintf('QUALITY METRICS'), 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.6, sprintf('Overall MSE: %.2e', overall_mse));
    text(0.1, 0.5, sprintf('SNR: %.1f dB', snr));
    text(0.1, 0.4, sprintf('Sparsity: %.1f%%', 100*sparsity_ratio));
    text(0.1, 0.3, sprintf('Spectral range: %.1f', range(mean_spectrum)));
    axis off;
    
    sgtitle('Reconstruction Quality Analysis', 'FontSize', 14, 'FontWeight', 'bold');
end

function rgb_composite = create_spectral_rgb_composite(hsi_cube, wavelengths)
    % Create RGB composite from hyperspectral cube
    [H, W, num_bands] = size(hsi_cube);
    
    % Find bands closest to standard RGB wavelengths
    target_red = 620;    % nm
    target_green = 540;  % nm
    target_blue = 460;   % nm
    
    [~, red_band] = min(abs(wavelengths - target_red));
    [~, green_band] = min(abs(wavelengths - target_green));
    [~, blue_band] = min(abs(wavelengths - target_blue));
    
    rgb_composite = zeros(H, W, 3);
    rgb_composite(:,:,1) = hsi_cube(:,:,red_band);
    rgb_composite(:,:,2) = hsi_cube(:,:,green_band);
    rgb_composite(:,:,3) = hsi_cube(:,:,blue_band);
    
    % Normalize for display
    for c = 1:3
        channel = rgb_composite(:,:,c);
        if max(channel(:)) > 0
            rgb_composite(:,:,c) = channel / max(channel(:));
        end
    end
end

% Keep all your existing FISTA functions exactly the same...
function x_recon = fista_lasso_diagonal(y, A, At, H, W, L, params)
    % FISTA algorithm for Lasso (your existing implementation)
    max_iter = params.max_iter;
    lambda = params.lambda;
    tol = params.tol;
    verbose = params.verbose;
    
    N = H * W * L;
    x = zeros(N, 1);
    z = x;
    t = 1;
    
    % Estimate Lipschitz constant
    if verbose
        fprintf('Estimating Lipschitz constant...\n');
    end
    Lipschitz = estimate_lipschitz_diagonal(A, At, N);
    
    % Precompute At(y) for efficiency
    Aty = At(y);
    
    if verbose
        fprintf('Lipschitz constant: %.6f\n', Lipschitz);
        fprintf('Running FISTA iterations:\n');
    end
    
    best_x = x;
    best_obj = inf;
    
    for iter = 1:max_iter
        x_prev = x;
        
        % Compute gradient: A'(A(z) - y) = A'(A(z)) - A'(y)
        Az = A(z);
        gradient = At(Az) - Aty;
        
        % Proximal gradient update
        x = soft_threshold(z - (1/Lipschitz) * gradient, lambda/Lipschitz);
        
        % FISTA acceleration
        t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
        z = x + ((t - 1) / t_next) * (x - x_prev);
        t = t_next;
        
        % Compute objective function
        if verbose && (mod(iter, 20) == 0 || iter == 1)
            residual = A(x) - y;
            data_fidelity = 0.5 * norm(residual)^2;
            regularization = lambda * norm(x, 1);
            objective = data_fidelity + regularization;
            
            rel_change = norm(x - x_prev) / (norm(x_prev) + eps);
            
            fprintf('Iter %4d: Obj = %.6e, DataFit = %.6e, Reg = %.6e, RelChange = %.6e\n', ...
                    iter, objective, data_fidelity, regularization, rel_change);
            
            % Track best solution
            if objective < best_obj
                best_obj = objective;
                best_x = x;
            end
        end
        
        % Check convergence
        if iter > 10
            rel_change = norm(x - x_prev) / (norm(x_prev) + eps);
            if rel_change < tol
                if verbose
                    fprintf('Converged at iteration %d with relative change %.6e\n', iter, rel_change);
                end
                x = best_x;
                break;
            end
        end
    end
    
    if iter == max_iter && verbose
        fprintf('Reached maximum iterations (%d)\n', max_iter);
        x = best_x;
    end
    
    % Reshape to hyperspectral cube and ensure non-negativity
    x_recon = reshape(x, [H, W, L]);
    x_recon = max(0, x_recon);
end

function y = soft_threshold(x, threshold)
    y = sign(x) .* max(abs(x) - threshold, 0);
end

function L = estimate_lipschitz_diagonal(A, At, N)
    x = randn(N, 1);
    x = x / norm(x);
    
    max_iter = 50;
    L_est = 0;
    
    for i = 1:max_iter
        y = A(x);
        x = At(y);
        L_new = norm(x);
        x = x / L_new;
        
        if abs(L_new - L_est) < 1e-6 * L_est
            break;
        end
        L_est = L_new;
    end
    
    L = L_est;
end