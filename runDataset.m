% runDataset
%   This file executes the mainSeg function across the entire Oasis brains
%   dataset, averages the scores, and prints the result.
%   Author: Tanya Nair
%   Last modified: Aug 13, 2016


%%
NCOMPONENTS = 3;
MAXITER_EM = 10;
MAXITER_ICM = 10;
STOPPERCENT = .5;

segType = {'otsu', 'kmeans'}; segIdx = 1;
BETA = [2 .66];     % weights for neighbourhood potentials
ALPHA = [0.5 .999 1.32];    % weights for unary potentials

%% Run one disc
oasisDir = ('C:\Users\tnair_000\Documents\medimage\OASIS\alldiscs\');
discDir = dir('C:\Users\tnair_000\Documents\medimage\OASIS\alldiscs\*O*');
score = cell(1,length(discDir));
for i=1:length(discDir)
    % brain masked, atlas registred, bias field corrected images
    imDir = [oasisDir discDir(i).name '\PROCESSED\MPRAGE\T88_111\'];
    hdr = dir([imDir '*_masked_gfc.hdr']);
    fprintf('%s\n', hdr.name); % sanity check
    hdrInfo = analyze75info([imDir hdr.name]);
    I = double(analyze75read( hdrInfo ));
    brainMask = I; brainMask(I==0)=1;

    % "ground truth" segmentations
    imDir = [oasisDir discDir(i).name '\FSL_SEG\'];
    hdr = dir([imDir '*_masked_gfc_fseg.hdr']);
    hdrInfo = analyze75info([imDir hdr.name]);
    Igt = double(analyze75read( hdrInfo ));
    IMDIMS = size(I);
    
    [I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType{segIdx}, brainMask);
    labels=I_initSeg;
    
    for iter=1:2
        [pxgn, labels, energy] = runICM( I, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA );
        [model, logli_img] = maximization_step( I, pxgn, model, NCOMPONENTS);
        I_finalSeg = reshape(labels,IMDIMS);
    end
    score{i} = scoreSeg(Igt, I_finalSeg, NCOMPONENTS, IMDIMS);
end
meanScore = mean(reshape(cell2mat(score),[4,3,length(discDir)]),3);

