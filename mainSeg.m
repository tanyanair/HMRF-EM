%% mainSeg
%  This file is the main executable that reads in a specified  image and
%  then performs HRMF-EM segmentation as described in Zhang et al. (2004)

%   Author: Tanya Nair
%   Last Modified: March 23, 2016

%% Constants
% close all; clear all;
segType = {'otsu', 'kmeans'}; segIdx = 1;

%% Load brain
imDir = 'C:\Users\Lakshmi\OneDrive\Documents\medimage\oasis_cross-sectional_disc1\disc1\OAS1_0001_MR1\PROCESSED\MPRAGE\T88_111';
hdr = [imDir '\OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc.hdr'];
hdrInfo = analyze75info(hdr);
I = double(analyze75read( hdrInfo ));
brainMask = I; brainMask(I==0)=1;

imDir = 'C:\Users\Lakshmi\OneDrive\Documents\medimage\oasis_cross-sectional_disc1\disc1\OAS1_0001_MR1\FSL_SEG\';
hdr = [imDir 'OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc_fseg.hdr'];
hdrInfo = analyze75info(hdr);
Igt = double(analyze75read( hdrInfo ));

betascale = .055;
% BETA = [0 .9 1.05 1];     % weights for neighbourhood potentials
% ALPHA = [0 1.1 1 1];    % weights for unary potentials
BETA = [1 1 1];     % weights for neighbourhood potentials
ALPHA = [1 1 1];    % weights for unary potentials
NCOMPONENTS = 3;
MAXITER_EM = 5;
MAXITER_ICM = 5;
STOPPERCENT = .5;
%%
IMDIMS = size(I);
% the model contains the mu & sig
[I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType{segIdx}, brainMask);

[I_finalSeg, model, energy, emIter] = runHMRF( I, I_initSeg, model, brainMask, NCOMPONENTS, ...
                                MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA, STOPPERCENT);

% I_finalSeg(brainMask==1)=1;
%% display results
figure;
subplot(131); imagesc( Igt(:,:,100) ); title( 'ground truth' );
subplot(132); imagesc( I(:,:,100) ); title( 'initial seg' );
subplot(133); imagesc( I_finalSeg(:,:,100) ); title( 'final seg' );
% title( sprintf( 'final seg, %i EM, %i ICM, %.3f beta ', MAXITER_EM, MAXITER_ICM, betascale )); 
           
figure;
subplot(121); plot(energy); title( 'Energy vs. EM Iteration' ); hold on; plot( energy, 'ro');
subplot(122); imagesc( I_finalSeg(:,:,100) ); title( 'final seg' );
% title( sprintf( 'final seg, %i EM, %i ICM, %.3f beta ', emIter, MAXITER_ICM, betascale ));

