clear
clc

RccTH = 0.95;
RsTH = 250;
RansacTH = 1e-9;
ransacRounds = 2000;
lrCorrTH = 5;
searchRadius = 100;
W = 3;

%% Read in images
image1Orig = imread('Cones_im2.jpg');
image1Gray = rgb2gray(image1Orig);
image2Orig = imread('Cones_im6.jpg');
image2Gray = rgb2gray(image2Orig);

%% Harris detector
Rs1 = harrisDetector(image1Gray, 100);
Rs2 = harrisDetector(image2Gray, 100);

%% NCC
image1 = double(image1Gray);
image2 = double(image2Gray);

filter = ones(7,7);
image1squared = image1 .* image1;
image2squared = image2 .* image2;

image1squared = imfilter(image1squared, filter, 'replicate');
image2squared = imfilter(image2squared, filter, 'replicate');

image1squared = image1squared.^0.5;
image2squared = image2squared.^0.5;

image1 = image1 ./image1squared;
image1(isnan(image1)) = 0;

image2 = image2 ./image2squared;
image2(isnan(image2)) = 0;

matches = struct();
count = 0;
for i = 4:(size(image1,1)-3)
    i
    for j = 4:(size(image1,2)-3)
        if Rs1(i,j) >= RsTH
            for k = 4:(size(image2,1)-3)
                for l = 4:(size(image2,2)-3)
                    if Rs2(k,l) >= RsTH
                        value = sum(sum(image1((i-3):(i+3),(j-3):(j+3)) .* image2((k-3):(k+3),(l-3):(l+3))));
                        if value >= RccTH
                            count = count + 1;
                            matches(count).row1 = i;
                            matches(count).col1 = j;
                            matches(count).row2 = k;
                            matches(count).col2 = l;
                            matches(count).val = value;
                        end
                    end
                end
            end
        end
    end
end

%% Compute Fundamental Matrix using RANSAC

%% NORMALIZE POINTS HERE??
% mean_x_1 = mean([matches.row1]);
% mean_y_1 = mean([matches.col1]);
% var1 = var([matches.row1 matches.row2], 0, 'all');
% T_1 = [var1^-1 0 -mean_x_1; 0 var1^-1 -mean_y_1; 0 0 1];
% points = [matches.row1; matches.col1; ones(size([matches.row1]))]
% for i = 1:size([matches.row1],2)
%     points(:,i) = points(:,i)' * T_1;
%     matches(i).row1 = points(1,i);
%     matches(i).col1 = points(2,i);
% end
%
% mean_x_2 = mean([matches.row2]);
% mean_y_2 = mean([matches.col2]);
% var2 = var([matches.row2 matches.row2], 0, 'all');
% T_2 = [var2^-1 0 -mean_x_2; 0 var2^-1 -mean_y_2; 0 0 1];
% points = [matches.row2; matches.col2; ones(size([matches.row2]))]
% for i = 1:size([matches.row2],2)
%     points(:,i) = points(:,i)' * T_2;
%     matches(i).row2 = points(1,i);
%     matches(i).col2 = points(2,i);
% end

%%
mean_x_1 = (size(image1,2)+1)/2;
mean_y_1 = (size(image1,1)+1)/2;
var_x_1 = var(1:size(image1,2));
var_y_1 = var(1:size(image1,1));
T_1 = [var_x_1^-1 0 -mean_x_1; 0 var_y_1^-1 -mean_y_1; 0 0 1];
points1 = [matches.col1; matches.row1; ones(size([matches.row1]))]
for i = 1:size([matches.row1],2)
    points1(:,i) = T_1 * points1(:,i);
    matches(i).col1 = points1(1,i);
    matches(i).row1 = points1(2,i);
end

mean_x_2 = (size(image2,2)+1)/2;
mean_y_2 = (size(image2,1)+1)/2;
var_x_2 = var(1:size(image1,2));
var_y_2 = var(1:size(image1,1));
T_2 = [var_x_2^-1 0 -mean_x_2; 0 var_y_2^-1 -mean_y_2; 0 0 1];
points2 = [matches.col2; matches.row2; ones(size([matches.row2]))]
for i = 1:size([matches.row2],2)
    points2(:,i) = T_2 * points2(:,i);
    matches(i).col2 = points2(1,i);
    matches(i).row2 = points2(2,i);
end
%%
usedPoints1 = struct('row',{},'col',{});
usedPoints2 = struct('row',{},'col',{});
pairs = struct();
vals = [matches.val];
totalPairings = 0;
for count = 1:size(matches,2)
    count
    [val, index] = max(vals);
    y1 = matches(index).row1;
    x1 = matches(index).col1;
    y2 = matches(index).row2;
    x2 = matches(index).col2;
    
    vals(index) = -1;
    
    used = false;
    for i = 1:size(usedPoints1,2)
        if(usedPoints1(i).row == y1 && usedPoints1(i).col == x1)
            used = true;
        end
    end
    for i = 1:size(usedPoints2,2)
        if(usedPoints2(i).row == y2 && usedPoints2(i).col == x2)
            used = true;
        end
    end
    
    if used
        continue
    end
    
    totalPairings = totalPairings + 1;
    
    usedPoints1(totalPairings).row = y1;
    usedPoints1(totalPairings).col = x1;
    usedPoints2(totalPairings).row = y2;
    usedPoints2(totalPairings).col = x2;
    
    pairs(totalPairings).row1 = y1;
    pairs(totalPairings).col1 = x1;
    pairs(totalPairings).row2 = y2;
    pairs(totalPairings).col2 = x2;
    pairs(totalPairings).val = val;
end

bestInlierCount = -1;
bestH = eye(3);
for num = 1:ransacRounds
    
    pointIndices = randi([1 size(pairs,2)], 1, 9);
    A = zeros(8,9);
    for i = 1:8
        A(i,:) = [pairs(pointIndices(i)).col1 * pairs(pointIndices(i)).col2,...
            pairs(pointIndices(i)).col1 * pairs(pointIndices(i)).row1, ...
            pairs(pointIndices(i)).col1,...
            pairs(pointIndices(i)).row1 * pairs(pointIndices(i)).col2, ...
            pairs(pointIndices(i)).row1 * pairs(pointIndices(i)).row2, ...
            pairs(pointIndices(i)).row1, ...
            pairs(pointIndices(i)).col2, ...
            pairs(pointIndices(i)).row2, 1];
    end
    
    [~,D,V] = svd(A);
    F = V(:,8);
    [U_f, D_f, V_f] = svd(F);
    [~, index] = min(D_f);
    D_f(index) = 0;
    F = U_f * D_f * transpose(V_f);
    F = transpose(T_2) * reshape(F, [3 3]) * T_1;
    F = reshape(F, [3 3]);
    inliersCount = 0;
    inliers = struct();
    
    for i = 1:size(pairs,2)
        y1 = pairs(i).row1;
        x1 = pairs(i).col1;
        y2 = pairs(i).row2;
        x2 = pairs(i).col2;
        p1 = [x1 y1 1];
        p2 = [x2 y2 1]';
        
        if abs(p1*F*p2)< RansacTH
            inliersCount = inliersCount + 1;
            inliers(inliersCount).row1 = pairs(i).row1;
            inliers(inliersCount).col1 = pairs(i).col1;
            inliers(inliersCount).row2 = pairs(i).row2;
            inliers(inliersCount).col2 = pairs(i).col2;
        end
    end
    
    if inliersCount > bestInlierCount
        bestInlierCount = inliersCount;
        bestF = F;
        bestInliers = inliers;
    end
end

% Compute a better F by using the closest 8 points
A = zeros(bestInlierCount,9);
for i = 1:bestInlierCount
    A(i,:) = [bestInliers(i).col1 * bestInliers(i).col2,...
        bestInliers(i).col1 * bestInliers(i).row1, ...
        bestInliers(i).col1,...
        bestInliers(i).row1 * bestInliers(i).col2, ...
        bestInliers(i).row1 * bestInliers(i).row2, ...
        bestInliers(i).row1, ...
        bestInliers(i).col2, ...
        bestInliers(i).row2, 1];
end

[~,D,V] = svd(A);
F = V(:,8);
[U_f, D_f, V_f] = svd(F);
[~, index] = min(D_f);
D_f(index) = 0;
F = U_f * D_f * transpose(V_f);
F = reshape(F, [3 3]);
F = transpose(T_2) * F * T_1;
points1 = [bestInliers.col1; bestInliers.row1; ones(size([bestInliers.row1]))];
for i = 1:bestInlierCount
    points1(:,i) = T_1 \ points1(:,i);
    bestInliers(i).col1 = points1(1,i);
    bestInliers(i).row1 = points1(2,i);
end

points2 = [bestInliers.col2; bestInliers.row2; ones(size([bestInliers.row2]))];
for i = 1:bestInlierCount
    points2(:,i) = T_2 \ points2(:,i);
    bestInliers(i).col2 = points2(1,i);
    bestInliers(i).row2 = points2(2,i);
end

%% Compute Disparity Maps

corr_window = 50;
search_window = 10;
%similarity_measure =

%%

xDisparity = zeros(size(image1));
yDisparity = zeros(size(image1));
maxY2 = (size(image2,1)-W);
maxY1 = (size(image1,1)-W);
for i = (1+W):(size(image1,1)-W)
    i
    for j = (1+W):(size(image1,2)-W)
        bestVal = Inf;
        bestK = -1;
        bestL = -1;
        abc = F * [j i 1]';
        a = -abc(1)/abc(2);
        b = -abc(3)/abc(2);
        image1Vals = image1Orig((i-W):(i+W),(j-W):(j+W),:);
        for l = max(1+W,j-searchRadius):min((size(image2,2)-W),j+searchRadius)
            k = round(a*l + b);
            if k >= (W+1) && k <= maxY2
                value = sum(sum(sum((image1Vals - image2Orig((k-W):(k+W),(l-W):(l+W),:)).^2)));
                if value < bestVal
                    bestVal = value;
                    bestL = l;
                end
            end
        end
        xDisparityTemp = bestL-j;
        bestK = a*bestL + b;
        yDisparityTemp = bestK-i;
        
        xDisparity(i,j) = xDisparityTemp;
        yDisparity(i,j) = yDisparityTemp;
        
        bestValTemp = Inf;
        bestJ = -1;
        kTemp = round(bestK);
        lTemp = bestL;
        abc = [lTemp kTemp 1] * F;
        a = -abc(1)/abc(2);
        b = -abc(3)/abc(2);
        image2Vals = image2Orig((kTemp-W):(kTemp+W),(lTemp-W):(lTemp+W),:);
        for jTemp = max(W+1,lTemp-searchRadius):min((size(image1,2)-W),lTemp+searchRadius)
            iTemp = round(a*jTemp + b);
            if iTemp >= (W+1) && iTemp <= maxY1
                value = sum(sum(sum((image1Orig((iTemp-W):(iTemp+W),(jTemp-W):(jTemp+W),:) - image2Vals).^2)));
                if value < bestValTemp
                    bestValTemp = value;
                    bestJ = jTemp;
                end
            end
        end
        
        if bestJ > -1
            if abs(bestJ - j) < lrCorrTH
                xDisparity(i,j) = xDisparityTemp;
                yDisparity(i,j) = yDisparityTemp;
            end
        end
    end
end

%% Normalize to images
xDisp = xDisparity- min(xDisparity(:));%abs(xDisparity);
yDisp = yDisparity- min(yDisparity(:));%abs(yDisparity);
xDisparityGrayImage = xDisp ./ max(xDisp(:));
yDisparityGrayImage = yDisp ./ max(yDisp(:)); 
%xDisparityGrayImage = medfilt2(xDisparityGrayImage);
%yDisparityGrayImage = medfilt2(yDisparityGrayImage);

figure(1);
imshow(xDisparityGrayImage)
figure(2);
imshow(yDisparityGrayImage)