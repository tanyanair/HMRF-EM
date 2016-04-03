function energy = calcEnergy2( labels, IMDIMS)

imLabels = reshape( labels, IMDIMS );
neighbourhoodMask = [ 0 1 0; 1 0 1; 0 1 0];
% neighbourhoodMask = [ 1 1 1; 1 0 1; 1 1 1];

labelmask = cell(1,2);
energy = cell(1,2);
for i=1:2
   labelmask{i} = zeros(IMDIMS);
   labelmask{i}(find(imLabels==i)) = 1;
   energy{i} = conv2(imLabels.*labelmask{i},neighbourhoodMask, 'same');
   energy{i} = energy{i}(:);
end
energy = -sum(cell2mat(energy),2);

