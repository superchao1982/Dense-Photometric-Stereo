function recsurf = buildModel(refinedNorms, texture)

% 
%  Using normals of all the pixels to get the reconstructed surface
%

%% Get depth of all pixels

s = size(refinedNorms);
slant = zeros(s(1:2));
tilt = zeros(s(1:2));
for i = 1:s(1)
    for j = 1:s(2)
        % change from image coordinate
        n_ij = squeeze(refinedNorms(s(1)+1-i,j,:));       
        % calculate slant & tilt
        [slant(i,j), tilt(i,j)] = grad2slanttilt(-n_ij(1)/n_ij(3), -n_ij(2)/n_ij(3));
    end
end

recsurf = shapeletsurf(slant, tilt, 6, 1, 2);
%figure, needleplotst(slant,tilt,5,3);
%% build and show the 3D model

[x, y] = meshgrid(1:s(2), 1:s(1));
figure('Name','Reconstructed Model'), ...
    h = surf(x,y,recsurf,'FaceColor',[0.498 0.816 0.922],'FaceAlpha',0.95,'EdgeColor','none');

%flip texture image for mapping
imgTexFlip = zeros(size(texture));
for i = 1:3
    imgTexFlip(:,:,i) = flipud(texture(:,:,i));
end
set(h, 'CData', imgTexFlip, 'FaceColor', 'texturemap');

camlight left;
lighting phong;
axis vis3d;
daspect([max(daspect)*[1 1] 2.5]);
axis off;

end