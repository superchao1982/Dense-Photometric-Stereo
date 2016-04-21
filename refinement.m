function refinedNorms = refinement(norms, lambda, sigma)
% 
%   description:  
%   Minimizing the energy function in MRF formulation by graph cut:
%   energy function E(N) = E_data(N) + E_smoothness(N),
%   where the data term is sum_s||N_s - N_{alpha_2}||,
%   "measure per-pixel differences between the measured and estimated ratio images"
%   and smoothness term is lambda * sum_{t is S's neighbor} log(1 + ||N_{alpha_s} - N_{alpha_t}|| / (2 * sigma * sigma))
%

%% create labels with z>0 on subdivided icosahedron

TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 5); 
vCoor = TR.X;
labels = vCoor(vCoor(:, 3) > 0, :); % labels corresponding to different normal orientations 
L = size(labels,1);  %|L| = 5057, same as stated in the paper

%% create and initialize a graph cut optimization object

s = size(norms);
gco = GCO_Create(s(1) * s(2), L);
% initialize data cost
% original E_data calculated is 0.XXXX, or 1.XXXX. Need a int32 type cost, 
% so *10000 to keep the precision, the same for E_smooth later 
E_data = int32( pdist2(labels, reshape(norms, [], 3)) * 10000 );                       
GCO_SetDataCost(gco, E_data);
% initialize smooth cost
dist = pdist2(labels, labels);
E_smooth = int32( lambda * log(1 + dist / (2 * sigma * sigma)) * 10000 );
GCO_SetSmoothCost(gco, E_smooth);

%% construct neighborhood matrix

value = ones(2*s(1)*s(2)-s(1)-s(2), 1);
idx_i = zeros(2*s(1)*s(2)-s(1)-s(2), 1);
idx_j = zeros(2*s(1)*s(2)-s(1)-s(2), 1);
k = 0;
for j = 1:s(2)
    for i = 1:s(1)
        if i < s(1)
            k = k+1;
            idx_i(k) = (j-1) * s(1) + i;
            idx_j(k) = idx_i(k) + 1;
        end
        if j < s(2)
            k = k+1;
            idx_i(k) = (j-1) * s(1) + i;
            idx_j(k) = idx_i(k) + s(1);
        end
    end
end
neighborM = sparse(idx_i, idx_j, value, s(1)*s(2), s(1)*s(2));
GCO_SetNeighbors(gco, neighborM);

%% get optimal labels and refined normal

GCO_Expansion(gco);
optLabels = GCO_GetLabeling(gco);
% reshape norms order from 1d to 2d (37976 ==> 202*188)
refinedNormsUnidxed = labels(optLabels, :);
refinedNorms = zeros(s);
for j = 1:s(2)
    for i = 1:s(1)
        refinedNorms(i,j,:) = refinedNormsUnidxed((j-1)*s(1) + i, :);
    end
end

% show the refined normal with L = (-1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
figure('Name','Refined Noraml'), ...
    imshow((-1/sqrt(3) * refinedNorms(:,:,1) + 1/sqrt(3) * refinedNorms(:,:,2) + 1/sqrt(3) * refinedNorms(:,:,3)) / 1.1);
end