clear, clc;

addpath('lib/gco/matlab','lib/S2_Sampling_Suite/S2_Sampling_Toolbox',...
    'lib/nearestneighbour','lib/surfPeterKovesi');
%GCO_UnitTest; % Get GCO lib initialized
SrcPath = 'data/data09';
SrcType = '*.bmp';

%% Resampling
[Imgs, LightVecs] = resampling(SrcPath,SrcType);

%% Normal initialization
[InitalNorms, idxDeImg] = normInit(Imgs, LightVecs);

%% Normal refinement

lambda = 0.5;    % weight of smoothness term
sigma = 0.65;    % smaller smoother

refinedNorms = refinement(InitalNorms, lambda, sigma);

%% Surface Constructio 
recsurf = buildModel(refinedNorms, Imgs(:,:,:,idxDeImg)/255);
