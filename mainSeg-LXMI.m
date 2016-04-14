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
% imDir = 'C:\Users\tnair_000\OneDrive\Documents\medimage\oasis_cross-sectional_disc1\disc1\OAS1_0001_MR1\PROCESSED\MPRAGE\T88_111';
hdr = [imDir '\OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc.hdr'];
hdrInfo = analyze75info(hdr);
I = double(analyze75read( hdrInfo ));
brainMask = I; brainMask(I==0)=1;

imDir = 'C:\Users\Lakshmi\OneDrive\Documents\medimage\oasis_cross-sectional_disc1\disc1\OAS1_0001_MR1\FSL_SEG\';
hdr = [imDir 'OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc_fseg.hdr'];
hdrInfo = analyze75info(hdr);
Igt = double(analyze75read( hdrInfo ));

% BETA = [0 .9 1.05 1];     % weights for neighbourhood potentials
% ALPHA = [0 1.1 1 1];    % weights for unary potentials
% BETA = [1.0101 1.0101 1.01];     % weights for neighbourhood potentials
NCOMPONENTS = 3;
MAXITER_EM = 10;
MAXITER_ICM = 10;
STOPPERCENT = .5;
IMDIMS = size(I);

%%
BETA = [2 .66];     % weights for neighbourhood potentials
ALPHA = [.53 .999 1.32];    % weights for unary potentials
[I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType{segIdx}, brainMask);
[I_finalSeg, model, energy, emIter] = runHMRF( I, I_initSeg, model, brainMask, NCOMPONENTS, ...
                                 MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA, STOPPERCENT);

score = [[0 ALPHA]' [0 BETA 0]' scoreSeg(Igt, I_finalSeg, NCOMPONENTS, IMDIMS)];
disp(score);

%% display results
figure;
subplot(131); imagesc( Igt(:,:,81) ); title( 'ground truth' );
subplot(132); imagesc( I_initSeg(:,:,81) ); title( 'initial seg' );
subplot(133); imagesc( I_finalSeg(:,:,81) ); title( 'final seg' ); 

% figure;
% subplot(121); plot(energy); title( 'Energy vs. EM Iteration' ); hold on; plot( energy, 'ro');
% subplot(122); imagesc( I_finalSeg(:,:,100) ); title( 'final seg' );
% title( sprintf( 'final seg, %i EM, %i ICM, %.3f beta ', emIter, MAXITER_ICM, betascale ));

%%

BETA = [2 .66];     % weights for neighbourhood potentials
ALPHA = [1 .99 1.32];    % weights for unary potentials
% the model contains the mu & sig
[I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType{segIdx}, brainMask);
labels=I_initSeg;
%
[pxgn, labels, energy] = runICM( I, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA );
[model, logli_img] = maximization_step( I, pxgn, model, NCOMPONENTS);

I_finalSeg = reshape(labels,IMDIMS);
score = [[0 ALPHA]' [0 BETA 0]' scoreSeg(Igt, I_finalSeg, NCOMPONENTS, IMDIMS)];
disp(score);

% display results
figure;
subplot(131); imagesc( Igt(:,:,81) ); title( 'ground truth' );
subplot(132); imagesc( I_initSeg(:,:,81) ); title( 'init seg' );
subplot(133); imagesc( I_finalSeg(:,:,81) ); title( sprintf('%.3f %.2', ALPHA, BETA) ); 
%
% figure; imagesc( I(:,:,90)); title( 'Original T1 MR0001 Slice 90' ); colormap jet
% figure; imagesc( Igt(:,:,90)); title( 'Ground Truth' ); colormap jet
%%
% i=1;
% for a=1:0.01:1.1
%    for b=1:0.01:1.1
%        for c=1:0.01:1.1
%             BETA = [c b a];     % weights for neighbourhood potentials
%             [I_finalSeg, model, energy, emIter] = runHMRF( I, I_initSeg, model, brainMask, NCOMPONENTS, ...
%                                 MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA, STOPPERCENT);
%             score = [[0 ALPHA]' [0 BETA]' scoreSeg(Igt, I_finalSeg, NCOMPONENTS, IMDIMS)];
%             disp(score);
%             res{i} = score;
%             i=i+1;
%        end
%    end
% end

