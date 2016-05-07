%% mainSeg
%  This file is the main executable that reads in a specified  image and
%  then performs HRMF-EM segmentation as described in Zhang et al. (2001)
%
%   Author: Tanya Nair
%   Last Modified: May 7, 2016

close all; clear all;

%% Load brain
imDir = 'C:\Users\tnair_000\OneDrive\Documents\medimage\oasis_cross-sectional_disc1\disc1\OAS1_0001_MR1\PROCESSED\MPRAGE\T88_111';
hdr = [imDir '\OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc.hdr'];
hdrInfo = analyze75info(hdr);
I = double(analyze75read( hdrInfo ));
brainMask = I; brainMask(I==0)=1;

imDir = 'C:\Users\tnair_000\OneDrive\Documents\medimage\oasis_cross-sectional_disc1\disc1\OAS1_0001_MR1\FSL_SEG\';
hdr = [imDir 'OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc_fseg.hdr'];
hdrInfo = analyze75info(hdr);
Igt = double(analyze75read( hdrInfo ));
NCOMPONENTS = 3;
MAXITER_EM = 10;
MAXITER_ICM = 10;
IMDIMS = size(I);

%%
BETA = [2 .66];     % weights for neighbourhood potentials
segType = 'otsu';
% segType = 'kmeans';
% ALPHA = [0.48 1.1 1.32];    % weights for unary potentials
ALPHA = [0.5 .999 1.32];
[I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType, brainMask);
[I_finalSeg, model, energy, emIter] = runHMRF( I, I_initSeg, model, brainMask, NCOMPONENTS, ...
                                 MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA);

score = [[0 ALPHA]' [0 BETA 0]' scoreSeg(Igt, I_finalSeg, NCOMPONENTS, IMDIMS)];
disp(score);

%% display results
figure;
slice = 83;
subplot(131); imagesc( Igt(:,:,slice) ); title( 'Ground Truth' );
subplot(132); imagesc( I_initSeg(:,:,slice) ); title( sprintf('Initial Segmentation: %s', segType) );
subplot(133); imagesc( I_finalSeg(:,:,slice) ); title( 'Final Segmentation' ); 

% figure;
% subplot(121); plot(energy(1:end-1)); title( 'Energy vs. EM Iteration' ); hold on; plot( energy(1:end-1), 'ro');
% subplot(122); imagesc( I_finalSeg(:,:,100) ); title( 'final seg' );
