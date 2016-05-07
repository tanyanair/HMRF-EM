% runDataset
%   This file executes the mainSeg function across the entire Oasis brains
%   dataset, averages the scores, and prints the result.
%
% Author: Tanya Nair
% Last Modified: May 7, 2016

close all; clear all;
%%
NCOMPONENTS = 3;
MAXITER_EM = 10;
MAXITER_ICM = 10;

% segType = 'otsu';
segType = 'kmeans';
BETA = [2 .66];     % weights for neighbourhood potentials
ALPHA = [0.5 .999 1.32];    % weights for unary potentials

oasisDir = ('E:\Datasets\medimage\OASIS\alldiscs\');
discDir = dir('E:\Datasets\medimage\OASIS\alldiscs\*O*');
scorefinal = cell(1,length(discDir));
scoreinit = cell(1,length(discDir));
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
    
    %%
    [I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType, brainMask);
    [I_finalSeg, model, energy, emIter] = runHMRF( I, I_initSeg, model, brainMask, NCOMPONENTS, ...
                                     MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA);
    scorefinal{i} = scoreSeg(Igt, I_finalSeg, NCOMPONENTS, IMDIMS);
    scoreinit{i} = scoreSeg(Igt, I_initSeg, NCOMPONENTS, IMDIMS);
    disp(scorefinal{i});
    

end
meanScoreFinal = mean(reshape(cell2mat(scorefinal),[4,3,length(discDir)]),3);
meanScoreInit = mean(reshape(cell2mat(scoreinit),[4,3,length(discDir)]),3);

