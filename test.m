% imageOrig = imread('DanaHallWay2/DSC_0286.JPG');
% imageGray = rgb2gray(imageOrig);
% figure; 
% % Rs = harrisDetector(imageGray, 100);
% figure
%  image('CData', image2Orig,'XData',[1 450], 'YData', [-1 -375])
% for i = 1:size(Rs2,1)
%     for j = 1:size(Rs2,2)
%         if(Rs2(i,j)>=150)
%             hold on
%             plot(j,-i,'.','MarkerSize',40)
%         end
%     end
% end
%
shiftNum = size(image2Orig,1);
shift = [1 0 0; 0 1 0; 0 shiftNum 1];
shiftHP = affine2d(shift);
[image1Shift image1ShiftRef] = imwarp(image1Orig, shiftHP);
image2Ref = imref2d(size(image2Orig));

figure(25);
clf;
imshowpair(image2Orig, image2Ref, image1Shift, image1ShiftRef, 'blend','Scaling','joint');

hold on

for i = 1:size(pairs,2)
    plot([pairs(i).col1, pairs(i).col2], [pairs(i).row1 + shiftNum, pairs(i).row2], 'Linewidth', 2);
end

figure(35);
clf;
imshowpair(image2Orig, image2Ref, image1Shift, image1ShiftRef, 'blend','Scaling','joint');

hold on

for i = 1:size(bestInliers,2)
    plot([bestInliers(i).col1, bestInliers(i).col2], [bestInliers(i).row1 + shiftNum, bestInliers(i).row2], 'Linewidth', 2);
end

% 

figure(45);
imshow(image1Orig);
hold on;

[x1 y1] = ginput(1);
x2 = 1:size(image2Orig,2);
abc = F * [x1 y1 1]';
a = -abc(1)/abc(2);
b = -abc(3)/abc(2);
y2 = a * x2 + b;

plot(x1,y1,'.','MarkerSize',20)

figure(55);
imshow(image2Orig)
hold on;
plot(x2, y2, 'Linewidth', 2)