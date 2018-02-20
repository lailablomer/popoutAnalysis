%% popoutAnalysis
% Takes a video file and returns the position of the mouse for every frame.
% Also it scores the time (in number of frames) the mouse spents at each of
% the 10 stimuli. If the tail of the mouse is found, the orientation is
% calculated for that position. All output is summarized in positions =
% [xposition, yposition, Area(0 if the mouse is no area),
% TailFound(logical), xtailposition, ytailposition]. 
%
% Code is available at https://github.com/lailablomer/popoutAnalysis 

function [positions, scoring, popout, rest] = popoutAnalysis(filename, startframe, endframe, draw, findCircle, includedArea, stimuli)

% stimuli = the amound of stimuli in the arena (normally 10, including the popout)

clc;    
close all;
workspace;

%% Retrieve file
folder = fileparts(which(filename));
movieFullFileName = fullfile(folder, filename);

% Check to see that it exists.
if ~exist(movieFullFileName, 'file')
	sprintf('File not found:\n%s\nRedefine filename and try again!', movieFullFileName);
	return;
end

numberOfFramesWritten = 0;

%% Background
% Create background by avereging frames and calculate the average maximum
% RGB values. 

videoObject = VideoReader(movieFullFileName);

% make background
[double_bg, numberOfFrames, meanMaxRGB] = makeBackground(startframe, endframe, videoObject);

%% Define areas
% Create n areas to score the video for. This can either be done by hand
% (draw = 1), the function can find a circle itself (findCircle = 1), or
% the user can define the circle diameter. 
[areaX, areaY, xOuter, yOuter] = defineAreas(draw, findCircle, stimuli, includedArea);

%% variables 
blackThreshold = 0.20; % Treshold for amount of blackness to be treated as mouse
n = stimuli; % amount of stimuli

oldpos = 0;
positions = zeros((numberOfFrames-startframe) + 1, 6);
mouseArea = zeros(1,10);

%% Video analysis
% For every frame, the background is substracted. Then, the resulting image is tresholded to have the remainig shape which is assumed to be the mouse. From this, the position of the mouse is calculated.
for frame = startframe : endframe
    
    vidFrames = read(videoObject, frame);
    Frame = vidFrames(:,:,:,1);

    % mouse head position is towards stimulus
    inArea = 0;
    area = 0;
    xBorderArea = [];

    % Subtraction and tresholding of current frame
    B = imcomplement(Frame);
    B = double(B) - double_bg;
    
    %% Find mouse position
    % Make a binary image of the mouse and check the RGB values. If the maximum RGB values are smaller than this value it is better to
    % have a normalized treshold since the video is pretty blurry.
    if meanMaxRGB > 90
        mouse = B(:,:,1) > meanMaxRGB* blackThreshold;
    else
        mouse = B(:,:,1) > max(B(:))* blackThreshold;
    end
    
    % get mouse with tail position
    [xtail, ytail, tail, finalpos] = mousePosition(oldpos, mouse, frame);
    oldpos = finalpos;

    %% Scoring
    % Check if mouse is in a specified area, and score if true. Also, check
    % if the mouse orientation is toward the stimulus.
    
    for i = 1:n
        if inpolygon(finalpos(1), finalpos(2), areaX(:,i), areaY(:,i))
            mouseArea(1,i) = mouseArea(1,i) + 1;
            xBorderArea = xOuter{i};
            yBorderArea = yOuter{i};
            area = i;
            break
        end
    end

    % check if mouse body position is in between tail and the border with
    % the stimulus. 
    if tail && ~isempty(xBorderArea)
        xTemp = horzcat(xtail ,xBorderArea', xtail);
        yTemp = horzcat(ytail, yBorderArea', ytail);
        inArea = inpolygon(finalpos(1), finalpos(2), xTemp, yTemp);
    end

    % save position
    positions(frame-startframe+1,:) = [finalpos, inArea, area, xtail, ytail];

    % keep track of processed frames
    if mod(frame+1, 100) == 0
        logmsg(['Processed frame ' num2str(frame-startframe+1) ' of ' num2str(numberOfFrames-startframe)]);
    end

    numberOfFramesWritten = numberOfFramesWritten + 1;
end

%% Collect data
scoring = struct('popout', mouseArea(1,10), 'area1', mouseArea(1,1), 'area2', mouseArea(1,2), 'area3', mouseArea(1,3), 'area4', mouseArea(1,4),'area5', mouseArea(1,5), 'area6', mouseArea(1,6), 'area7', mouseArea(1,7), 'area8', mouseArea(1,8), 'area9', mouseArea(1,9));
popout = mouseArea(1,10);
rest = sum(mouseArea(1,1:9)) / 9;

% alert user when done
logmsg(['Done!  It processed ' num2str(numberOfFramesWritten), ' frames of ' movieFullFileName]);
logmsg('Circles and popout area are drawn over background in figure 1');

end
