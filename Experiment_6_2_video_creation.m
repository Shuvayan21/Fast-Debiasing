%% Create compact comparison video from 4 folders without white space

clc; clear; close all;

% --- Step 1: Define folder paths and parameters ---
folderNames = {
    'Original', ...
    'LASSO', ...
    'Debiased LASSO', ...
    'coded_snapshots_grayscale'
};

folderTitles = {
    'Original', ...
    'LASSO', ...
    'Debiased LASSO', ...
    'coded_snapshots_grayscale'
};

% Assuming 30 sets with 5 frames each
numSets = 30;
framesPerSet = 5;

% --- Step 2: Verify all folders exist ---
fprintf('Checking folder existence...\n');
for i = 1:length(folderNames)
    if ~exist(folderNames{i}, 'dir')
        error('Folder does not exist: %s', folderNames{i});
    else
        fprintf('✓ Folder found: %s\n', folderNames{i});
    end
end

% --- Step 3: Get file lists and verify file counts ---
fprintf('\nChecking file counts in each folder...\n');

% For first three folders: should have 30 sets × 5 frames = 150 files
expectedFilesFirstThree = numSets * framesPerSet;
fileLists = cell(1, 3);

for i = 1:3
    files = dir(fullfile(folderNames{i}, '*.png'));
    fileLists{i} = files;
    
    if length(files) ~= expectedFilesFirstThree
        warning('Folder %s has %d files, expected %d', ...
            folderNames{i}, length(files), expectedFilesFirstThree);
    else
        fprintf('✓ Folder %s: %d files\n', folderNames{i}, length(files));
    end
end

% For coded_snapshots_grayscale: should have 30 files
codedFiles = dir(fullfile(folderNames{4}, '*.png'));
if length(codedFiles) ~= numSets
    warning('Folder %s has %d files, expected %d', ...
        folderNames{4}, length(codedFiles), numSets);
else
    fprintf('✓ Folder %s: %d files\n', folderNames{4}, length(codedFiles));
end

% --- Step 4: Sort files in proper order ---
fprintf('\nSorting files...\n');

% Sort coded snapshots by set number
codedFilePaths = cell(numSets, 1);
for setIdx = 1:numSets
    filename = sprintf('coded_snapshot_set%02d_gray.png', setIdx);
    filepath = fullfile(folderNames{4}, filename);
    if exist(filepath, 'file')
        codedFilePaths{setIdx} = filepath;
    else
        error('Coded snapshot file not found: %s', filename);
    end
end

% Sort files for first three folders by set and frame
sortedFilePaths = cell(3, numSets, framesPerSet);
for folderIdx = 1:3
    for setIdx = 1:numSets
        for frameIdx = 1:framesPerSet
            filename = sprintf('set%02d_frame%02d.png', setIdx, frameIdx);
            filepath = fullfile(folderNames{folderIdx}, filename);
            if exist(filepath, 'file')
                sortedFilePaths{folderIdx, setIdx, frameIdx} = filepath;
            else
                error('File not found: %s', filename);
            end
        end
    end
    fprintf('✓ Sorted files for: %s\n', folderNames{folderIdx});
end

% --- Step 5: Read sample image to get dimensions ---
sampleImage = imread(sortedFilePaths{1, 1, 1});
[height, width, ~] = size(sampleImage);
fprintf('\nImage dimensions: %d x %d\n', height, width);

% --- Step 6: Create video writer ---
outputVideoName = 'compact_comparison_video.mp4';
videoWriter = VideoWriter(outputVideoName, 'MPEG-4');
videoWriter.FrameRate = 2; % 2 frames per second
videoWriter.Quality = 95; % High quality

fprintf('\nCreating video: %s\n', outputVideoName);
open(videoWriter);

% --- Step 7: Generate video frames without white space ---
fprintf('Generating compact video frames...\n');
totalFrames = numSets * framesPerSet;
frameCount = 0;

% Calculate dimensions for side-by-side layout
numCols = 4; % 4 methods
totalWidth = width * numCols;
totalHeight = height;

for setIdx = 1:numSets
    fprintf('Processing set %d/%d...\n', setIdx, numSets);
    
    % Get coded snapshot for this set
    codedImage = imread(codedFilePaths{setIdx});
    
    for frameIdx = 1:framesPerSet
        frameCount = frameCount + 1;
        
        % Read images from first three folders
        images = cell(1, 3);
        for folderIdx = 1:3
            images{folderIdx} = imread(sortedFilePaths{folderIdx, setIdx, frameIdx});
        end
        
        % Create combined frame without white space
        sampleCombined = zeros(totalHeight, totalWidth, 'uint8');
        
        % Place images side by side
        combinedFrame(:, 1:width) = images{1}; % Original
        combinedFrame(:, width+1:2*width) = images{2}; % FISTA
        combinedFrame(:, 2*width+1:3*width) = images{3}; % Deb
        combinedFrame(:, 3*width+1:4*width) = codedImage; % Coded Snapshot
        
        % Add text labels directly on the images
        labeledFrame = insertText(combinedFrame, [width*0.5-50, 20], folderTitles{1}, ...
            'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
        labeledFrame = insertText(labeledFrame, [width*1.5-50, 20], folderTitles{2}, ...
            'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
        labeledFrame = insertText(labeledFrame, [width*2.5-50, 20], folderTitles{3}, ...
            'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
        labeledFrame = insertText(labeledFrame, [width*3.5-50, 20], folderTitles{4}, ...
            'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
        
        % Add frame and set information
        infoText = sprintf('Set %02d - Frame %02d/%d', setIdx, frameIdx, framesPerSet);
        labeledFrame = insertText(labeledFrame, [totalWidth-200, totalHeight-40], infoText, ...
            'FontSize', 18, 'TextColor', 'white', 'BoxColor', 'red', 'BoxOpacity', 0.8, 'AnchorPoint', 'RightCenter');
        
        % Write frame to video
        writeVideo(videoWriter, labeledFrame);
        
        % Display progress
        if mod(frameCount, 5) == 0
            fprintf('  Progress: %d/%d frames processed\n', frameCount, totalFrames);
        end
    end
end

% --- Step 8: Close video writer ---
close(videoWriter);

% --- Step 9: Display completion message and video info ---
fprintf('\n🎉 Compact comparison video creation complete!\n');
fprintf('   Output file: %s\n', outputVideoName);
fprintf('   Total frames: %d\n', totalFrames);
fprintf('   Video duration: %.1f seconds\n', totalFrames / videoWriter.FrameRate);
fprintf('   Frame rate: %d fps\n', videoWriter.FrameRate);
fprintf('   Final resolution: %d x %d\n', totalWidth, totalHeight);
fprintf('   Layout: 4 images side-by-side without gaps\n');

% --- Step 10: Show sample frame from the video ---
fprintf('\nDisplaying sample frame...\n');
figure('Position', [100, 100, totalWidth/2, totalHeight/2]);

% Create a sample frame from set 1, frame 1
sampleImages = cell(1, 3);
for folderIdx = 1:3
    sampleImages{folderIdx} = imread(sortedFilePaths{folderIdx, 1, 1});
end
sampleCoded = imread(codedFilePaths{1});

sampleCombined = zeros(totalHeight, totalWidth, 'uint8');
sampleCombined(:, 1:width) = sampleImages{1};
sampleCombined(:, width+1:2*width) = sampleImages{2};
sampleCombined(:, 2*width+1:3*width) = sampleImages{3};
sampleCombined(:, 3*width+1:4*width) = sampleCoded;

% Add labels for the sample
sampleLabeled = insertText(sampleCombined, [width*0.5-50, 20], folderTitles{1}, ...
    'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
sampleLabeled = insertText(sampleLabeled, [width*1.5-50, 20], folderTitles{2}, ...
    'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
sampleLabeled = insertText(sampleLabeled, [width*2.5-50, 20], folderTitles{3}, ...
    'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');
sampleLabeled = insertText(sampleLabeled, [width*3.5-50, 20], folderTitles{4}, ...
    'FontSize', 20, 'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.7, 'AnchorPoint', 'Center');

sampleLabeled = insertText(sampleLabeled, [totalWidth-200, totalHeight-40], 'Set 01 - Frame 01/05', ...
    'FontSize', 18, 'TextColor', 'white', 'BoxColor', 'red', 'BoxOpacity', 0.8, 'AnchorPoint', 'RightCenter');

imshow(sampleLabeled);
title('Sample Frame from Compact Comparison Video', 'FontSize', 14, 'FontWeight', 'bold');

% --- Step 11: Verify video file was created ---
if exist(outputVideoName, 'file')
    videoInfo = dir(outputVideoName);
    fileSizeMB = videoInfo.bytes / (1024^2);
    fprintf('\nVideo file created successfully!\n');
    fprintf('   File size: %.2f MB\n', fileSizeMB);
    
    % Offer to play the video
    choice = questdlg('Would you like to play the created video?', ...
        'Video Complete', 'Yes', 'No', 'Yes');
    if strcmp(choice, 'Yes')
        try
            winopen(outputVideoName);
            fprintf('   Opening video...\n');
        catch
            fprintf('   Could not automatically open the video.\n');
            fprintf('   Please manually open: %s\n', outputVideoName);
        end
    end
else
    warning('Video file was not created successfully.');
end

fprintf('\nCompact comparison video process completed successfully!\n');