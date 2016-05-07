%% getInitSeg returns an initial segmentation of IMG into NCOMPONENTS 
%   using the SEGTYPE specified.  The SEGTYPE must be either otsu or
%   kmeans. The function also returns the mean, variance, and inverse
%   frequency for each component.
%
% Author: Tanya Nair
% Last Modified: May 7, 2016

function [seg, model] = getInitSeg( img, ncomponents, segType, brainMask)

switch segType
    case 'kmeans'
        % Approximate initial means are used to ensure consistent labeling
        % of the segmentation. These values are rounded off from the means
        % obtained through Otsu thresholding.
        seg = kmeans(img(:), ncomponents, 'start', [300; 700; 1000] );
        seg = reshape(seg,size(img));
    case 'otsu'
        seg = cell(1,size(img,3));
        for i=1:size(img,3)
            seg{i} = otsu(img(:,:,i), ncomponents);
        end
        seg = reshape(cell2mat(seg),size(img));
end
seg(brainMask==1)=0;

mu = zeros( ncomponents,1 );
sigma = zeros( ncomponents,1 );
alpha = zeros( ncomponents, 1);
for i=1:ncomponents
    img_components = img(seg==i);
    mu(i) = mean(img_components);
    sigma(i) = var(img_components);
    alpha(i) = 1/ncomponents;
end
model.mu = mu;
model.sig = sigma;
model.alpha = alpha;
