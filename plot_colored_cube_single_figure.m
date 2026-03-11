function plot_colored_cube_single_figure(data, mymax, wavelengths, indices)
% PLOT_COLORED_CUBE_SINGLE_FIGURE Plots all spectral slices in a single figure window
% Input:
%   data = datacube to plot (height x width x num_spectrums)
%   mymax = scaling factor to increase brightness of each slice
%   wavelengths = wavelengths of each spectral slice in nm
%   indices = which spectral indices to plot

    data = double(data);
    data = mymax * data / max(data(:));
    
    num_plots = length(indices);
    
    % Calculate subplot grid dimensions
    cols = 7;
    rows = 4;
    
    % Create single figure window
    figure('Name', 'Hyperspectral Data Cube - All Wavelengths', ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 800]);
    
    for i = 1:num_plots
        r = indices(i);
        currentwl = wavelengths(r);
        
        % Create subplot
        subplot(rows, cols, i);
        
        % Create colormap for this wavelength
        cmM = gray * kron(ones(3, 1), spectrumRGB(wavelengths(r)));
        cmM = cmM / max(max(cmM));
        
        % Display the image
        subimage(squeeze(data(:,:,r)), colormap(cmM)); 
        axis off;
        axis image;
        
        % Add title with wavelength information
        title(sprintf('%.1f nm\n(Band %d)', currentwl, r), ...
              'FontSize', 9, 'FontWeight', 'bold');
    end
    
    % Add overall title
    sgtitle('Hyperspectral Image - All Wavelengths', 'FontSize', 14, 'FontWeight', 'bold');
end