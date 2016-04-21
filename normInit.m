function [InitalNorm, deImgIdx] = normInit(Imgs, LightVecs)

%   Note: Imgs dimension: [size(img) #images],  Imgs(:,:,:,i) is the ith image;
%
%   description:
%
%   1. Selection of the denominator image (least affected by shdows and hightlights): 
%      See paper section 4.3 
%      "For each resampled image_i, count the number of pixels whose intensity rank
%      satisfies rank_i > L, where L > median. Let k_L^i be the total number of pixels in image_i 
%      satisfying this condition, r_L^i be the mean rank among the pixels in image_i that satisfies
%      this condition. The denominator image is defined to be one with 
%      1) maximum k_L and 2) r_L lower than some thresh H. 
%      Currently, we set L and H to be the 70th and 90th percentiles respectively."
%    
%
%   2. Initial local normal estimation by ratio images
%      Every pixel in an ratio image is expressed by:
%               I_1     N * L_1
%              ----- = ---------
%               I_2     N * L_2
%
%      can be transform into form of A_x*N_x + A_y*N_y + A_z*N_z = 0
%      where A(k) = I_1 * L_2(k) - I_2*L_1(k)
%      take denominator image as reference, construct k-1 equations,
%      ==> Solve AX = 0 by SVD, explicitly enforces ||N|| = 1
%   

%% Denominator Image

ImgNum = size(Imgs,4);
s = [size(Imgs,1) size(Imgs,2)];
L = 0.7 * ImgNum;
H = 0.9 * ImgNum;

% get pixel intensities
grayImgs = zeros([s ImgNum]);
for i = 1:ImgNum
    grayImgs(:,:,i) = 0.2989 * Imgs(:,:,1,i) + 0.5870 * Imgs(:,:,2,i) + 0.1140 * Imgs(:,:,3,i);
end
% rank pixel intensities, loop over all pixels
rankPixel = zeros([s ImgNum]);
rankIJ = zeros(ImgNum, 1);
for i = 1:s(1)
    for j = 1:s(2)
        [~, ind] = sort(grayImgs(i,j,:)); % get sorted index
        rankIJ(ind) = 1:ImgNum; % reverse sort
        rankPixel(i,j,:) = rankIJ;
    end
end

% thresholding on L and H to get denominator image index
info = zeros(ImgNum, 2);
for i = 1:ImgNum
    % count pixels with rank > L in each image
    highRankInd = (rankPixel(:,:,i) >= L); 
    info(i,1) = sum(highRankInd(:));
    % calculate mean rank among the pixels that satisfy above consition
    highRankPx = rankPixel(:,:,i) .* highRankInd;
    info(i,2) = sum(highRankPx(:)) / info(i,1);
end
[~, deImgIdx] = max(info(:,1) .* (info(:,2) < H));

figure('Name','denominator image'), imshow(Imgs(:,:,:,deImgIdx)/255);
%figure('Name','denominator image intensity'), imshow(grayImgs(:,:,deImgIdx)/255);

%% Initial Normal Estimation

%light vectors for the denominator image and the others
lv_deImg = LightVecs(deImgIdx, :);
lv_rest = [LightVecs(1:(deImgIdx-1),:); LightVecs((deImgIdx+1):end,:)];

% calculate the initial normal
InitalNorm = zeros([s 3]);
for i = 1:s(1)
    for j = 1:s(2)
        pxIJ = squeeze(grayImgs(i,j,:));
        I_deImg = pxIJ(deImgIdx);
        I_rest = [pxIJ(1:(deImgIdx-1)); pxIJ((deImgIdx+1):end)];
        A = I_rest * lv_deImg - I_deImg * lv_rest;  
        [~,~,V] = svd(A,0);
        % reverse those with negative z
        InitalNorm(i,j,:) = sign(V(3,3))*V(:,3); 
    end
end
% show the recovered normal with L = (-1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
figure('Name','Initial Noraml'), ...
    imshow((-1/sqrt(3) * InitalNorm(:,:,1) + 1/sqrt(3) * InitalNorm(:,:,2) + 1/sqrt(3) * InitalNorm(:,:,3)) / 1.1);

end