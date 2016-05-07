%% calcEnergy3D calculates the energy of the image labels for NCOMPONENTS
%   using the BETA scaling, a first-order 3D neighbourhood (6 voxels), and
%   a delta-dirac as the energy potential function.
%
% Author: Tanya Nair
% Last Modified: May 7, 2016

function energy = calcEnergy3D( labels, brainMask, IMDIMS, NCOMPONENTS, BETA)

imLabels = reshape( labels, IMDIMS );
% 3D first order neighbourhood mask
neighbourhoodMask = zeros(3,3,3);
neighbourhoodMask(:,:,1) = BETA(2)*[0 0 0; 0 1 0; 0 0 0];
neighbourhoodMask(:,:,2) = BETA(1)*[0 1 0; 1 0 1; 0 1 0];
neighbourhoodMask(:,:,3) = BETA(2)*[0 0 0; 0 1 0; 0 0 0];

labelmask = cell(1,NCOMPONENTS);
energy = cell(1,NCOMPONENTS);
for i=1:NCOMPONENTS
   labelmask{i} = zeros(IMDIMS);
   labelmask{i}(imLabels==i) = 1;
   energy{i} = convn(imLabels.*labelmask{i}/i, neighbourhoodMask, 'same' );
   energy{i} = -energy{i}(:);
   energy{i}(brainMask==1)=0;
end
