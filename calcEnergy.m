function energy = calcEnergy( labels, IMDIMS)

imLabels = reshape( labels, IMDIMS );
imLabels_pad = padarray(imLabels,[1,1]);

N = cell(1,4); % this is the neighbourhood
N{1} = imLabels_pad(1:end-2,2:end-1); % north
N{2} = imLabels_pad(3:end,2:end-1); % south
N{3} = imLabels_pad(2:end-1,3:end); % east
N{4} = imLabels_pad(2:end-1,1:end-2); % west

for i=1:length(N)
    N{i} = N{i}(:);
end
N = cell2mat(N);

potentials = calcPotentials(labels, N);
energy = sum(potentials,2); % the sum of the potentials


