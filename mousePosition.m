%% mousePosition
% This function takes the binary image with the mouse and computes the
% position of the mouse, and, if possible the tail position. 

function [xtail, ytail, tail, finalpos] = mousePosition(oldpos, mouse, frame)
% Set variables
minAreaSize = 200; % Minimal area size for region that is tracked as mouse
beginFound = 0;
tailNotFound = 0;
tailWidth = 12; 
tailToMiddle = 70;

%% Find mouse position
% If the position is not found, either the old position is taken as the
% current one. Worste case scenario there is no shape, and the position is
% set to be 0.
pos = regionprops(mouse, 'Centroid', 'Area');
if isempty(pos)
    if oldpos ~= 0
        finalpos = oldpos;
    else
        logmsg(['!!Cannot find position in frame: ', frame,'!!']);
        finalpos = [0 0];
        xtail = 0;
        ytail = 0;
        tail = 0;

        return;
    end
else
    % Check whether the areasize of the found shape is larger than the
    % minimum
    maxAreaInd = find([pos.Area] == max([pos.Area]));
    maxAreaInd = maxAreaInd(1);
    nearMaxInd = ([pos.Area] > pos(maxAreaInd).Area * 0.2 & [pos.Area] > minAreaSize);
    
    % take previous position as mouse position and alert user
    if pos(maxAreaInd).Area <= minAreaSize
        logmsg(['Could not find position in frame: ', frame, ', take previous position as new position']);
        finalpos = oldpos;
        
    % else, take centre of found shape as mouse position
    else
        posCentroids = [pos(nearMaxInd).Centroid];
        finalpos = [mean(posCentroids(1:2:end)), mean(posCentroids(2:2:end))];
    end
end

%% Mouse boundaries
% Get mouse boundaries and make sure there is only one shape (the mouse) in the image
boundary = bwboundaries(mouse);
A = cellfun('size', boundary, 1);
[~, ind] = max(A);
mouseBoundary = boundary{ind};

% design new binary image with 1 shape, the mouse. Also make new
% mouseBoundaries
[M, N] = size(mouse);
mouseBinary = poly2mask(mouseBoundary(:,2), mouseBoundary(:,1), M, N);
mouseBinary = bwmorph(mouseBinary, 'bridge', Inf);

mouseBinary = bwmorph(mouseBinary, 'thicken');
boundary = bwboundaries(mouseBinary);
A = cellfun('size', boundary, 1);
[~, ind] = max(A);
mouseBoundary = boundary{ind};
[row, ~] = size(mouseBoundary);

%% Find tail
% find farthest geodesic point from mouse position and check if it is far
% enough from mouse centre.
D = bwdistgeodesic(mouseBinary, floor(finalpos(1)), floor(finalpos(2)), 'quasi-euclidean');
posTails = D(:);

% Define distance between tail and body and check length

% % dist = tailToMiddle - 1;
% % while dist < tailToMiddle
% %     [num,indD] = max(posTails);
% %     if isnan(num) || (num == 0)
% %         tailNotFound = 1;
% %         break
% %     end
% %     [ytail,xtail] = ind2sub(size(D),indD);
% %     posTails(indD) = 0;
% %     dist = pdist([finalpos; [xtail, ytail]]);    
% % end

[~,indD] = max(posTails);
[ytail,xtail] = ind2sub(size(D),indD);
dist = pdist([finalpos; [xtail, ytail]]); 
if dist < tailToMiddle
    tailNotFound = 1;
end
    
% Compute distance between found tail and every point of mouseBoundary get the coordinate with the minimum distance
if ~tailNotFound
    np = length(mouseBoundary(:,1));
    Pp = [ytail xtail];

    % matrix of distances between all points and all vertices
    dpv(:,:) = hypot((repmat(mouseBoundary(:,1)', [np 1])-repmat(Pp(:,1), [1 1])),...
                     (repmat(mouseBoundary(:,2)', [np 1])-repmat(Pp(:,2), [1 1])));

    % Find the vector of minimum distances to vertices. 
    [~, index] = min(abs(dpv),[],2);
    ind = index(1);
    firstInd = ind;
    left = ind;
    right = ind;
    halfway = round(mod(ind + (row / 2), row));
end

%% Beginning of tail
% Look for beginning of tail by following the mouse boundary untill
% distance between the two sides become larger than tailWidth.
while ~beginFound && ~tailNotFound
    left = mod((left + row - 2), row) + 1;
    right = mod(right, row) + 1;
    dist = pdist([mouseBoundary(left,:); mouseBoundary(right,:)]);
    
    if dist > tailWidth
        beginFound = 1;
        ind = right;
    elseif (left == halfway) || (right == halfway)
        tailNotFound = 1;
    end
end

% get tail coordinates
if tailNotFound
    tail = 0;
    xtail = 0;
    ytail = 0;
else
    tail = 1;
    
%% Separate tail and body
% Check if the end of the mouseBoundary is passed or not. Based on this,
% devide the mouse in the binary image in the body and the tail

    passEnd = 0;
    if ind > firstInd
        dif = ind - firstInd;
        tailB = mod(firstInd - dif, row);
        tailE = ind;
        if tailB == 0
            tailB = 1;
        end
    else
        dif = firstInd - ind;
        tailB = ind;
        tailE = mod(firstInd + dif, row);
        if tailE == 0
            tailE = row;
        end
    end
    
    if tailB > tailE
        passEnd = 1;
    end
    
    % separate tail and body
    if passEnd
        fulltail = vertcat(mouseBoundary(tailB:end, :),mouseBoundary(1:tailE, :));
        rest = mouseBoundary(tailE:tailB, :);
    else
        fulltail = mouseBoundary(tailB:tailE, :);
        rest = vertcat(mouseBoundary(tailE:end, :), mouseBoundary(1:tailB, :));
    end

%% Find new tail and mouse positions

    % new tail position
    xtail = round(fulltail(1,2) + fulltail(end,2)) / 2;
    ytail = round(fulltail(1,1) + fulltail(end,1)) / 2;
    
    % new mouse binary mask
    tailBinary = poly2mask(fulltail(:,2), fulltail(:,1), M, N);
    mouseNew = poly2mask(rest(:,2), rest(:,1), M, N);
    
    % if there are multiple shapes in the binary image, make a new image
    % with only the largest shape 
    temp = bwconncomp(mouseNew);
    if temp.NumObjects > 1
        numPixels = cellfun(@numel,temp.PixelIdxList);
        [~ ,idx] = max(numPixels);
        temp2 = zeros(size(mouse));
        temp2(temp.PixelIdxList{idx}) = 1;
        mouseNew = logical(temp2);
    end
    
    % only if the body is bigger than the tail, take the new mouse position
    % from the centre of the body. 
    if (sum(tailBinary(:)) < sum(mouseNew(:)))
        posNew = regionprops(mouseNew, 'Centroid');
        finalpos = posNew.Centroid;
    end
end
end
