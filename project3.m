clear
clc

RccTH = 0.9;
RsTH = 50;
RansacTH = 2;
ransacRounds = 300;

%% Read in images
image1Orig = imread('DanaHallWay2/DSC_0285.JPG');
image1Gray = rgb2gray(image1Orig);
image2Orig = imread('DanaHallWay2/DSC_0286.JPG');
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
    
    pointIndices = randi([1 size(pairs,2)], 1, 4);
    
    point1_0 = [pairs(pointIndices(1)).col1 pairs(pointIndices(1)).row1];
    point2_0 = [pairs(pointIndices(2)).col1 pairs(pointIndices(2)).row1];
    point3_0 = [pairs(pointIndices(3)).col1 pairs(pointIndices(3)).row1];
    point4_0 = [pairs(pointIndices(4)).col1 pairs(pointIndices(4)).row1];
    
    point1_1 = [pairs(pointIndices(1)).col2 pairs(pointIndices(1)).row2];
    point2_1 = [pairs(pointIndices(2)).col2 pairs(pointIndices(2)).row2];
    point3_1 = [pairs(pointIndices(3)).col2 pairs(pointIndices(3)).row2];
    point4_1 = [pairs(pointIndices(4)).col2 pairs(pointIndices(4)).row2];
    
    A = zeros(8,8);
    A(1,:) = [point1_0 1 0 0 0 (-point1_0*point1_1(1))];
    A(2,:) = [0 0 0 point1_0 1 (-point1_0*point1_1(2))];
    A(3,:) = [point2_0 1 0 0 0 (-point2_0*point2_1(1))];
    A(4,:) = [0 0 0 point2_0 1 (-point2_0*point2_1(2))];
    A(5,:) = [point3_0 1 0 0 0 (-point3_0*point3_1(1))];
    A(6,:) = [0 0 0 point3_0 1 (-point3_0*point3_1(2))];
    A(7,:) = [point4_0 1 0 0 0 (-point4_0*point4_1(1))];
    A(8,:) = [0 0 0 point4_0 1 (-point4_0*point4_1(2))];
    
    b = [point1_1 point2_1 point3_1 point4_1]';
    
    hFlat = A\b;
    
    h = [hFlat(1:3)'; hFlat(4:6)'; hFlat(7:8)' 1];
    
    inliersCount = 0;
    inliers = struct();
    
    for i = 1:size(pairs,2)
        y1 = pairs(i).row1;
        x1 = pairs(i).col1;
        y2 = pairs(i).row2;
        x2 = pairs(i).col2;
        p1 = [x1 y1 1]';
        p2 = h * p1;
        p2 = p2./p2(3);
        p2 = round(p2);
        if abs(y2 - p2(2))<RansacTH && abs(x2 - p2(1))<RansacTH
            inliersCount = inliersCount + 1;
            inliers(inliersCount).row1 = pairs(i).row1;
            inliers(inliersCount).col1 = pairs(i).col1;
            inliers(inliersCount).row2 = pairs(i).row2;
            inliers(inliersCount).col2 = pairs(i).col2;
        end
    end
    
    if inliersCount > bestInlierCount
        bestInlierCount = inliersCount;
        bestH = h;
        bestInliers = inliers;
    end
end

A = zeros(bestInlierCount*2,8);
b = zeros(bestInlierCount*2,1);
for i = 1:bestInlierCount
    point1_0 = [bestInliers(i).col1 bestInliers(i).row1];
    point1_1 = [bestInliers(i).col2 bestInliers(i).row2];
    A(2*i-1,:) = [point1_0 1 0 0 0 (-point1_0*point1_1(1))];
    A(2*i,:) = [0 0 0 point1_0 1 (-point1_0*point1_1(2))];
    b(2*i-1,1) = point1_1(1);
    b(2*i,1) = point1_1(2);
    hFlat = A\b;
    
    totalH = [hFlat(1:3)'; hFlat(4:6)'; hFlat(7:8)' 1];
end
%%
%test = [1 0 0; 0 1 0; -75 0 1];
%testHP = affine2d(test);
%[testTF testRef] = imwarp(image1Orig, testHP);

image2Ref = imref2d(size(image2Orig));

totalHP = projective2d(totalH');
[tfImage tfImageRef] = imwarp(image1Orig,totalHP);
figure;
clf;
imshowpair(image2Orig,image2Ref, tfImage, tfImageRef, 'blend','Scaling','joint')
%imshowpair(image2Orig, image2Ref, testTF, testRef,'blend','Scaling','joint')