%% defineAreas
% This function shapes the areas that will be used to count the position of
% the mouse. The areas can be defined automatically by drawing a circle
% (draw = 0) or the user can draw the areas by hand on the background
% figure (draw = 1).

function [areaX, areaY, xOuter, yOuter] = defineAreas(draw, findCircle, n, includedArea)
% Setting parameters
outerDiameter = 500;

%% Draw circle by hand
% If draw is true the user should draw every area (n in total) by hand
% on the figure, these areas will be counted when defining the mouse
% position. The last area should be the popout area.
if draw
    title('Draw the popout LAST');
        
    for i = 1:n        
        H = imfreehand;
        binaryImage = H.createMask();
        A = bwboundaries(binaryImage);
        Axy = A{1}; % Get n by 2 array of x,y coordinates
        areaX(i) = Axy(:, 2); % Columns
        areaY(i) = Axy(:, 1); % Rows
        hold on;
    end

%% Define circle automatically
% If the circle should be calculated, first try to automatically detect a
% circle in the background, otherwise the diameter has to be set by hand.
else 
    noCircle = 0;
    
    % automatically find circle
    if findCircle
        [center, radius] = imfindcircles(backgroundImage, [(outerDiameter/2) ((outerDiameter/2)+100)], 'Sensitivity', 0.9);
        
        % if no circle is found, alert the user and ask to define the
        % diameter. 
        if isempty(center)
            logmsg(['Could not find outer circle, re-select diameter!']); 
            noCircle = 1;
        end
    end
  
    %% Select diameter by hand
    % If no circle can be found, the user has to set the diameter by hand. From this line the centre and radius will be defined.
    if ~noCircle        
        title('SELECT DIAMETER!');
        line1 = imdistline(gca);
        pause;
        
        % get the radius
        api1 = iptgetapi(line1);
        radius = round(api1.getDistance() / 2);
        
        % define middle of the circle
        pos = api1.getPosition;
        center = [mean(pos(:,1)) mean(pos(:,2))];
    end
    
    logmsg(['Circle found!']);
    
    %% Draw circle
    % define middle and radia of the outer and inner circle
    xCenter = center(1);
    yCenter = center(2);
    radiusOuter = radius;
    radiusInner = radiusOuter - (includedArea / 2);
    
    %% Select popout position
    % define position of popout by clicking on the figure
    title('Select popout position');
    [Popx, Popy] = getpts();
    hold on;
    logmsg('Popout position is defined');
    
    %% Define areas
    
    % Outer circle dynamics
    theta = 0 : 0.03 : 2*pi;
    xOuter = radiusOuter * cos(theta) + xCenter;
    yOuter = radiusOuter * sin(theta) + yCenter;
    lengthOuter = length(xOuter);

    % Inner circle dynamics
    if includedArea == 0
        xInner(1,1:lengthOuter) = xCenter;
        yInner(1,1:lengthOuter) = yCenter;
    else
        xInner = radiusInner * cos(theta) + xCenter;
        yInner = radiusInner * sin(theta) + yCenter;
    end
    lengthInner = length(xInner);
    
    % find popout position on circle 
    xDif = xOuter - Popx;
    yDif = yOuter - Popy;
    [~, ind] = min(abs(xDif - yDif));
    
    % If the positions are too far apart
    if ~(yOuter(ind) > Popy - 100) && (yOuter(ind) < Popy + 100)
        ind = mod(ind + (lengthOuter / 2), lengthOuter);
    end
    
    % Re-define start and end to be at the popout 
    ind = mod(round(ind + (lengthOuter/(n*2))), lengthOuter);
    xOuter = [xOuter(ind:end), xOuter(1:ind-1)];
    yOuter = [yOuter(ind:end), yOuter(1:ind-1)];
    xInner = [xInner(ind:end), xInner(1:ind-1)];
    yInner = [yInner(ind:end), yInner(1:ind-1)];
    
    % Plot circles on the figure
    plot(xOuter, yOuter);
    hold on;
    plot(xInner, yInner);
    hold on;
    
    %% Devide circles in n areas
    
    % Reshape values in n cells
    xOuter = num2cell(reshape(round(xOuter), lengthOuter/n, n),1);
    yOuter = num2cell(reshape(round(yOuter), lengthOuter/n, n),1);

    xInner = num2cell(reshape(round(xInner), lengthInner/n, n),1);
    yInner = num2cell(reshape(round(yInner), lengthInner/n, n),1);
%%
    [len, ~] = size(xOuter{1});
    
    areaX = zeros(2*len + 1, n);
    areaY = zeros(2*len + 1, n);

    % Define x and y values for n parts of the ring
    for i = 1:n
        areaX(:,i) = [xOuter{i}; flip(xInner{i}); xOuter{i}(1)];
        areaY(:,i) = [yOuter{i}; flip(yInner{i}); yOuter{i}(1)];
    end

    % plot popout-position on figure
    plot(areaX(:,10),areaY(:,10));
    hold on;
end
            
end