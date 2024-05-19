% Initialize webcam
cam = webcam;
% cam.Resolution='640x480';

% Initialize frame counter and total processing time
frameCount = 0;
totalProcessingTime = 0;

% Parameters
sigma = 1;
kernel_size = 5;
weak_pixel = 75;
strong_pixel = 255;
lowThresholdRatio = 0.02;
highThresholdRatio = 0.1;

% Create figure for displaying the original and edge-detected video
figureHandle = figure;
% Set the figure to full screen
set(figureHandle, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

while frameCount < 50
    % Capture a frame from the webcam
    input_image = snapshot(cam);
    frameCount = frameCount + 1;

    % Start timing the processing
    tic;

    % Convert the frame to grayscale
    grayImage = rgb2gray(input_image);
    
    % Convert the grayscale image to double precision
    img = imgaussfilt(grayImage, sigma);
    img = double(img);
    
    % Sobel Operator masks
    Kx = [-1 0 1; -2 0 2; -1 0 1];
    Ky = [1 2 1; 0 0 0; -1 -2 -1];

    % Gradient approximations
    Ix = imfilter(img, double(Kx), 'conv', 'replicate');
    Iy = imfilter(img, double(Ky), 'conv', 'replicate');
    
%     Ix = imgaussfilt(Ix, sigma);
%     Iy = imgaussfilt(Iy, sigma);

    % Gradient magnitude and direction
    G_mag = sqrt(Ix.^2 + Iy.^2);
    G_mag = G_mag / max(G_mag(:)) * 255; % Normalization
    angle = atan2(Iy, Ix)*180/pi;
    
    [height, width] = size(G_mag);

    % Non-maximum suppression
    output = zeros(height, width);
    angle(angle < 0) = angle(angle < 0) + 180;

    for i = 2:height-1
        for j = 2:width-1
            q = 255;
            r = 255;

            % Angle 0
            if (0 <= angle(i, j) && angle(i, j) < 22.5) || (157.5 <= angle(i, j) && angle(i, j) <= 180)
                q = G_mag(i, j + 1);
                r = G_mag(i, j - 1);
            % Angle 45
            elseif (22.5 <= angle(i, j) && angle(i, j) < 67.5)
                q = G_mag(i + 1, j - 1);
                r = G_mag(i - 1, j + 1);
            % Angle 90
            elseif (67.5 <= angle(i, j) && angle(i, j) < 112.5)
                q = G_mag(i + 1, j);
                r = G_mag(i - 1, j);
            % Angle 135
            elseif (112.5 <= angle(i, j) && angle(i, j) < 157.5)
                q = G_mag(i - 1, j - 1);
                r = G_mag(i + 1, j + 1);
            end

            if (G_mag(i, j) >= q) && (G_mag(i, j) >= r)
                output(i, j) = G_mag(i, j);
            else
                output(i, j) = 0;
            end
        end
    end

    % Double threshold
    highThreshold = max(output(:)) * highThresholdRatio;
    lowThreshold = highThreshold * lowThresholdRatio;

    res = zeros(height, width);

    strong = strong_pixel;
    weak = weak_pixel;

    strong_indices = find(output >= highThreshold);
    weak_indices = find((output < highThreshold) & (output >= lowThreshold));

    res(strong_indices) = strong;
    res(weak_indices) = weak;

    % Hysteresis
    for i = 2:height-1
        for j = 2:width-1
            if (res(i, j) == weak)
                if any(any(res(i-1:i+1, j-1:j+1) == strong))
                    res(i, j) = strong;
                else
                    res(i, j) = 0;
                end
            end
        end
    end

    % Stop measuring processing time
    processingTime = toc;
    totalProcessingTime = totalProcessingTime + processingTime;

    %Display original video and edge-detected video in one figure
    subplot(1, 2, 1);
    imshow(input_image);
    title('Original Video');

    subplot(1, 2, 2);
    imshow(res);
    title('Edge-detected Video');

    %Allow display to refresh
    drawnow;
end

clear cam;

% Calculate average processing time per frame
averageProcessingTime = totalProcessingTime / frameCount;

% Display performance metrics
fprintf('Total Processing Time: %.4f seconds\n', totalProcessingTime);
fprintf('Average Processing Time per Frame: %.4f seconds\n', averageProcessingTime);
fprintf('Frame Rate: %.2f FPS\n', 1 / averageProcessingTime);

disp('Processing complete!');
