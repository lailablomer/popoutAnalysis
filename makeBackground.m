%% makeBackground
% Function makes background from videoObject from the given start- until the
% endframe, and it calculates the illuminance through the maximum RGB
% values. 

function [double_bg, numberOfFrames, meanMaxRGB] = makeBackground(startframe, endframe, videoObject)

% if endframe is zero the whole videoObject from the startframe will be
% taken
if endframe == 0
    numberOfFrames = videoObject.NumberOfFrames;
else
    numberOfFrames = endframe;
end

figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

%% Average frames
% The background is complemented so black shapes become white and can be
% substracted from each other. 
bgframes = (startframe:90:numberOfFrames);
firstdone = 0;

% loop through frames and sum
for i = bgframes
    vidFrames = read(videoObject, i);
    Frame=vidFrames(:,:,:,1);
    if ~firstdone
        bgsum = double(Frame);
        firstdone = 1;
    else
        bgsum = bgsum + double(Frame);
    end
end

% devide by number of frames
bg = bgsum/length(bgframes);
double_bg = double(imcomplement(uint8(bg)));

% show image
imshow(uint8(double_bg));
hold on;
logmsg(['Background is made']);

%% Calculate max RGB values
% This is important since different videos might have different illuminance
% levels
testFrames = bgframes;
maxRGBs = zeros(1,length(testFrames));
j = 1;

% loop through frames and take max RGB value of every frame
for i = testFrames
    vidFrames = read(videoObject, i);
    Frame = vidFrames(:,:,:,1);
    B = imcomplement(Frame);
    B = double(B) - double_bg;
    maxRGBs(j) = max(B(:));
    j = j + 1;
end

% calculate the mean
meanMaxRGB = mean(maxRGBs);
logmsg(['Maximum RGB values are calculated']);
end