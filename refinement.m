function refinedNorms = refinement(norms, lambda, sigma)
% 
%   description:  
%   Minimizing the energy function in MRF formulation by graph cut:
%   energy function E(N) = E_data(N) + E_smoothness(N),
%   where the data term is sum_s||N_s - N_{alpha_2}||,
%   "measure per-pixel differences between the measured and estimated ratio images"
%   and smoothness term is lambda * sum_{t is S's neighbor} log(1 + ||N_{alpha_s} - N_{alpha_t}|| / (2 * sigma * sigma))
%
% ######## STEPS ABOUT using GCO library: my reference is ########
% ######## http://blog.csdn.net/on2way/article/details/46742765 ########
%
%% create labels with z>0 on subdivided icosahedron

TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 5); 
vCoor = TR.X;
labels = vCoor(vCoor(:, 3) > 0, :); % labels corresponding to different normal orientations 
L = size(labels,1);  %|L| = 5057, same as stated in the paper

%% create and initialize a graph cut optimization object

size_n = size(norms);
handle = GCO_Create(size_n(1) * size_n(2), L);
% 1. set data cost
% original E_data calculated is 0.XXXX, or 1.XXXX. Need a int32 type cost, 
% so *10000 to keep the precision, the same for E_smooth later 
E_data = int32(pdist2(labels, reshape(norms, [], 3)) * 10000);                       
GCO_SetDataCost(handle, E_data);
% 2. set neighbors
% construct a sparse 0-1 neighborhood matrix
value = ones(2*size_n(1)*size_n(2)-size_n(1)-size_n(2), 1);
idx_i = zeros(2*size_n(1)*size_n(2)-size_n(1)-size_n(2), 1);
idx_j = zeros(2*size_n(1)*size_n(2)-size_n(1)-size_n(2), 1);
k = 0;
% organize the indices
for i = 1:size_n(1)
    for j = 1:size_n(2)
        if i < size_n(1)
            k = k+1;
            idx_i(k) = i + (j-1) * size_n(1);
            idx_j(k) = idx_i(k) + 1;
        end
        if j < size_n(2)
            k = k+1;
            idx_i(k) = i + (j-1) * size_n(1);
            idx_j(k) = idx_i(k) + size_n(1);
        end
    end
end
weights = sparse(idx_i, idx_j, value, size_n(1)*size_n(2), size_n(1)*size_n(2));
GCO_SetNeighbors(handle, weights);

% 3. set smooth cost
E_smooth = int32(lambda * log(1 + pdist2(labels, labels) / (2 * sigma * sigma)) * 10000);
GCO_SetSmoothCost(handle, E_smooth);

%% get optimal labels and refined normal

GCO_Expansion(handle);
bestLabels = GCO_GetLabeling(handle);
% reshape norms order from 1d to 2d (37976 ==> 202*188)
refinedNormsUnidxed = labels(bestLabels, :);
refinedNorms = zeros(size_n);
for i = 1:size_n(1)
    for j = 1:size_n(2)
        refinedNorms(i,j,:) = refinedNormsUnidxed(i + (j-1)*size_n(1), :);
    end
end

% show the refined normal with L = (-1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
figure('Name','Refined Noraml'), ...
    imshow((-1/sqrt(3) * refinedNorms(:,:,1) + 1/sqrt(3) * refinedNorms(:,:,2) + 1/sqrt(3) * refinedNorms(:,:,3)) / 1.1);
end