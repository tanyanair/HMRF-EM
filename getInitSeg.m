function [seg, model] = getInitSeg( img, ncomponents, segType, brainMask)

% should do otsu thresholding here, for now just use kmeans
switch segType
    case 'kmeans'
        seg = kmeans(img',ncomponents);
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
