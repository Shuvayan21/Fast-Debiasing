%% Extract 30 sets of 5 frames from a video with top 100 rows removed

clc; clear;

% --- Step 1: Load video ---
videoFile = 'waterfall.mp4';   % Replace with your video filename
v = VideoReader(videoFile);
totalFrames = floor(v.Duration * v.FrameRate);

% --- Step 2: Define parameters ---
numSets = 30;       % number of sets
framesPerSet = 5;   % frames per set
rowsToRemove = 100; % remove top 100 rows

% --- Step 3: Compute frame indices ---
% Use first available frames for 30 sets
maxStartFrame = totalFrames - framesPerSet + 1;
setStarts = 1:framesPerSet:min(framesPerSet * numSets, maxStartFrame);
numSets = min(numSets, length(setStarts)); % Adjust if not enough frames

% --- Step 4: Create output directory ---
outputDir = 'extracted_frames_30';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- Step 5: Initialize storage ---
% Get frame dimensions after cropping
v.CurrentTime = 0;
sampleFrame = readFrame(v);
[originalHeight, originalWidth, ~] = size(sampleFrame);
newHeight = originalHeight - rowsToRemove;
frameSize = [newHeight, originalWidth];

videoSets = zeros([frameSize, 3, framesPerSet, numSets], 'uint8'); 

% --- Step 6: Loop through sets and extract frames ---
for s = 1:numSets
    startFrame = setStarts(s);
    
    for f = 1:framesPerSet
        frameIdx = startFrame + f - 1;
        
        % Set current time and read frame
        v.CurrentTime = (frameIdx - 1) / v.FrameRate;
        frame = readFrame(v);
        
        % Remove top 100 rows
        frame_cropped = frame(rowsToRemove+1:end, :, :);
        
        % Store in variable
        videoSets(:, :, :, f, s) = frame_cropped;
        
        % Save as individual file
        filename = sprintf('%s/set%02d_frame%02d.png', outputDir, s, f);
        imwrite(frame_cropped, filename);
        
        fprintf('Saved: %s\n', filename);
    end
    
    fprintf('✅ Completed set %d/%d\n', s, numSets);
end

% --- Step 7: Save the variable as well ---
save('extracted_video_sets.mat', 'videoSets', 'setStarts', 'numSets', 'framesPerSet', '-v7.3');

% --- Step 8: Display summary ---
fprintf('\n🎉 Extraction Complete!\n');
fprintf('   Total sets extracted: %d\n', numSets);
fprintf('   Frames per set: %d\n', framesPerSet);
fprintf('   Total frames: %d\n', numSets * framesPerSet);
fprintf('   Original frame size: %d x %d\n', originalHeight, originalWidth);
fprintf('   Cropped frame size: %d x %d\n', newHeight, originalWidth);
fprintf('   Files saved in: %s/\n', outputDir);
fprintf('   Variable saved as: extracted_video_sets.mat\n');

% --- Optional: Display sample frame from first set ---
figure('Position', [100, 100, 800, 400]);
subplot(1,2,1);
imshow(frame); title('Original Frame (first frame)');
subplot(1,2,2);
imshow(videoSets(:,:,:,1,1)); title('Cropped Frame (top 100 rows removed)');


