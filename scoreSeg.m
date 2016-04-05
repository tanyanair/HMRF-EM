function scores = scoreSeg( Igt, Iseg, NCOMPONENTS, IMDIMS)

%% similarity metric: percentage of correct classifications (not including background)
ncorrect = zeros(1,NCOMPONENTS+1); 
dice = zeros(1,NCOMPONENTS+1);
sensitivity = zeros(1,NCOMPONENTS+1);
specificity = zeros(1,NCOMPONENTS+1);
for k=0:NCOMPONENTS
   labelmaskGT = zeros(IMDIMS);
   labelmaskSeg = zeros(IMDIMS);
   
   labelmaskGT(Igt==k) = 1;
   labelmaskSeg(Iseg==k) = 1;
   
   ncorrect(k+1) = nnz(labelmaskGT.*labelmaskSeg);
   
   % sensitivity = true positive rate
   sensitivity(k+1) = ncorrect(k+1) / nnz(labelmaskGT);

   % specificity = false negative rate
   t0 = double(~labelmaskGT);
   p0 = double(~labelmaskSeg);
   specificity(k+1) = nnz(p0.*t0) / nnz(t0);
   
   % dice score
   dice(k+1) = ncorrect(k+1)*2 / (nnz(labelmaskSeg) + nnz(labelmaskGT));   
end
%% print these results
scores = [dice' sensitivity' specificity'];
% disp(res);
