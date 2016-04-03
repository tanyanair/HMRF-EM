function energy = calcEnergy3D( labels, IMDIMS, NCOMPONENTS)

imLabels = reshape( labels, IMDIMS );
% 2D
% neighbourhoodMask = [ 0 1 0; 1 0 1; 0 1 0];

% 3D mask
neighbourhoodMask = zeros(3,3,3);
neighbourhoodMask(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
neighbourhoodMask(:,:,2) = [0 1 0; 1 0 1; 0 1 0];
neighbourhoodMask(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

%% TODO: CHECK IF THIS WORKS
labelmask = cell(1,NCOMPONENTS);
energy = cell(1,NCOMPONENTS);
for i=1:NCOMPONENTS
   labelmask{i} = zeros(IMDIMS);
   labelmask{i}(find(imLabels==i)) = 1;
%    energy{i} = conv2(imLabels.*labelmask{i},neighbourhoodMask, 'same');
   energy{i} = convn(imLabels.*labelmask{i}, neighbourhoodMask, 'same' );
   energy{i} = energy{i}(:);
end
energy = -sum(cell2mat(energy),2);

