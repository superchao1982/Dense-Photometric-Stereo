function [ResampledImg, L_o] = resampling(SrcPath, SrcType)

%  Do uniform resampling on images
%
%  See paper section 4.2:
%  "seek the nearest light di-rection L_o at one vertex in the 
%  subdivided icosahedron for each captured light direction L_i, 
%  and interpolate the image I_o at L_o by:
%
%                                                 L_o * L_i' * I_i(x,y)
%  I_o(x,y) = sum_{i | L_i's NN == L_o} * -------------------------------------
%                                          sum_{i | L_i's NN == L_o} L_o * L_i'
% 
% "

%% find nearest light direction L_o of input light direction L_i on the subdivided icosahedron

% read captured light direction                            
L_i = textread([SrcPath '/lightvec.txt']);
L_i = normr(L_i);
% create the subdivided icosahedron, get vertices             
TR = IcosahedronMesh;  
TR = SubdivideSphericalMesh(TR, 4);
vCoor = TR.X;
% find the index of nearest neighbour in icosahedron vertices
NNIndex = nearestneighbour(L_i', vCoor'); 
% get rid of those repeated vertices
[uniqueIndex, ~, revIndex] = unique(NNIndex);
uniqueNum = length(uniqueIndex); 
L_o = zeros(uniqueNum, 3);
for i = 1:length(uniqueIndex)
    L_o(i,:) = vCoor(uniqueIndex(i),:,:);
end

%% Interpotate image I_o at L_o

ImgsInfo = dir([SrcPath '/' SrcType]);
ImgNum = size(ImgsInfo, 1);
ResampledImg = zeros([size(double(imread([SrcPath '/' ImgsInfo(1).name]))) uniqueNum]);
Weight = zeros(uniqueNum, 1); % the denominator of the equation

% calculate the equation
% numerator
for i = 1:ImgNum
    cur_image = double(imread([SrcPath '/' ImgsInfo(i).name]));
    ri = revIndex(i);
    w = vCoor(NNIndex(i),:,:) * L_i(i,:,:)';
    ResampledImg(:,:,:,ri) = ResampledImg(:,:,:,ri) + w * cur_image;
    Weight(ri) = Weight(ri) + w;
end
% divide by the denominator
for i = 1:uniqueNum
    ResampledImg(:,:,:,i) = ResampledImg(:,:,:,i) / Weight(i);
end

end