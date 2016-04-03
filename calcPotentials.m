function potential = calcPotentials(imLabels, N)

imLabels=imLabels(:);
% px_neighbours = [imLabels N];
potential_mtx = zeros(numel(imLabels),size(N,2));
for i=1:numel(imLabels)
    for j=1:size(N,2)
        potential_mtx(i,j) = dirac(imLabels(i) - N(i,j));
    end
end
potential_mtx(potential_mtx(:,:)==Inf) = 1;
potential = sum(potential_mtx,2);