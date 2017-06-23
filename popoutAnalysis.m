%% mousePosition 
% Takes a video file and returns the position of the mouse for every frame.
% also it scores the time (in number of frames) the mouse spents at each of
% the 10 stimuli.

function [positions, scoring, popout, rest, files] = popoutAnalysis(filename, startframe, endframe, draw, findCircle, includedArea)
%% 
% filename = file to be analysed
% startframe, endframe = both integers
% draw = true or false: specifies whether you want to draw the ROI's yourself. 
% excludedArea = integer, specifies the part of a circle excluded from analysis. 
%   If 0, the whole outer circle is analysed.

% tic;

% clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
% clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

%% variables 
blackThreshold = 0.20; % Treshold for amount of blackness to be treated as mouse %was 0.30
n = 10; % amount of stimuli

%% Retrieve file
folder = fileparts(which(filename));
movieFullFileName = fullfile(folder, filename);
% Check to see that it exists.
if ~exist(movieFullFileName, 'file')
	sprintf('File not found:\n%s\nRedefine filename and try again!', movieFullFileName);
	return;
end

numberOfFramesWritten = 0;

% try
    %% background
    videoObject = VideoReader(movieFullFileName);
	
    % make background
    [double_bg, numberOfFrames, meanMaxRGB] = makeBackground(startframe, endframe, videoObject);
    
    %% define boundaries around stimuli
    % draw boundary on first frame
    [areaX, areaY, xOuter, yOuter] = defineAreas(draw, findCircle, n, includedArea);
    
%% 
    oldpos = 0;
    positions = zeros((numberOfFrames-startframe) + 1, 6);
    mouseArea = zeros(1,10);
    
    count = 1;
    %% actual analysis part
    for frame = startframe : endframe		
        % For every frame, the background is substracted. Then, the
        % resulting image is tresholded to have the remainig shape which is assumed
        % to be the mouse. From this, the position of the mouse is calculated.
    %% 
        vidFrames = read(videoObject, frame);
        Frame = vidFrames(:,:,:,1);
        
        % mouse head position is towards stimulus
        inArea = 0;
        area = 0;
        xBorderArea = [];
        
        % Subtraction and tresholding of current frame
        B = imcomplement(Frame);
        B = double(B) - double_bg;
        
%         if mod(frame, 100) == 0
%             figure
%             imshow(B);
%         end
        
        
        %% find mouse position
        % If the maximum RGB values are smaller than this value it is better to
        % have a normalized treshold since the video is pretty blurry.
        
        if meanMaxRGB > 90
            mouse = B(:,:,1) > meanMaxRGB* blackThreshold;
        else
            mouse = B(:,:,1) > max(B(:))* blackThreshold;
        end
        
%         if mod(frame, 100) == 0
%             figure
%             imshow(mouse);
            filename = strcat('frame', num2str(frame));
            files{count} = struct('image', mouse, 'frame', filename);
            count = count+1;
%             imsave(mouse);
%         end
        
        % get mouse with tail position
        [xtail, ytail, tail, finalpos] = mousePosition(oldpos, mouse, frame);
        oldpos = finalpos;
        
%         if mod(frame, 250) == 0
%            figure;
%            imshow(mouse);
%            hold on;
%            if tail
%                plot(xtail,ytail,'*r');
%                hold on;
%            end
%                
%            plot(finalpos(1), finalpos(2),'*g');
%            hold on;
%         end
        
        %% track frames spend at stimulus
        % check if mouse is in a specified area
        for i = 1:n
            if inpolygon(finalpos(1), finalpos(2), areaX(:,i), areaY(:,i))
                mouseArea(1,i) = mouseArea(1,i) + 1;
                xBorderArea = xOuter{i};
                yBorderArea = yOuter{i};
                area = i;
                break
            end
        end
        
        % check if mouse body position is in between tail and border
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
        
%         if numberOfFramesWritten == 50
%             return
%         end
    end
    
%% collect data
    scoring = struct('popout', mouseArea(1,10), 'area1', mouseArea(1,1), 'area2', mouseArea(1,2), 'area3', mouseArea(1,3), 'area4', mouseArea(1,4),'area5', mouseArea(1,5), 'area6', mouseArea(1,6), 'area7', mouseArea(1,7), 'area8', mouseArea(1,8), 'area9', mouseArea(1,9));
    popout = mouseArea(1,10);
    rest = sum(mouseArea(1,1:9)) / 9;
	
	% alert user when done
	logmsg(['Done!  It processed ' num2str(numberOfFramesWritten), ' frames of ' movieFullFileName]);
	logmsg('Circles and popout area are drawn over background in figure 1');
    
% catch ME
% 	% Some error happened if you get here.
% 	strErrorMessage = sprintf('Error extracting movie frames from:\n\n%s\n\nError: %s\n\n)', movieFullFileName, ME.message);
% 	uiwait(msgbox(strErrorMessage));
% end

% toc;

end
