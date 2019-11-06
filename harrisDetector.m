function [Rs] = harrisDetector(image, nonMaxThreshold)

% Calculate gradient masks
[Ix,Iy] = imgradientxy(image, 'prewitt');

% calculate gradient products
IxIx = Ix.*Ix;
IyIy = Iy.*Iy;
IxIy = Ix.*Iy;


% Sum in 3x3 window
filter = ones(3,3);
IxIx = imfilter(IxIx, filter, 'replicate');
IxIy = imfilter(IxIy, filter, 'replicate');
IyIy = imfilter(IyIy, filter, 'replicate');

% Compute C matrix and R matrix
k = .04;
Rs = zeros(size(Ix));
for i = 1:size(Ix,1)
    for j = 1:size(Ix,2)
        C = [IxIx(i,j) IxIy(i,j); IxIy(i,j) IyIy(i,j)];
        eigs = eig(double(C));
        R = eigs(1)*eigs(2) - (k * (eigs(1) + eigs(2))^2);
        Rs(i,j) = R;
    end
end

%% Nonmax suppression
 
[Imag, Idir] = imgradient(image, 'prewitt');
Imagcells = mat2cell(Imag, [75 75 75 75 75], [90 90 90 90 90]);
cells = mat2cell(Rs,[75 75 75 75 75], [90 90 90 90 90]); 
for i = 1:5
    for j = 1:5
        vec = reshape(cells{i,j},1,[]); 
        [~, index] = sort(vec, 'descend'); 
        vec = zeros(75, 90); 
        vec(index(1:30)) = Imagcells{i,j}(index(1:30));
         for k = 30:-1:1
            if min(abs(index(1:k-1) - index(k))) < 10 | index(k) < nonMaxThreshold 
                vec(index(k)) = 0; 
            end 
         end 
        cells{i,j} = vec; 
    end 
end 
Rs = cell2mat(cells); 
