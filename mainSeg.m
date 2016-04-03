%% mainSeg
%  This file is the main executable that reads in a specified  image and
%  then performs HRMF-EM segmentation as described in Zhang et al. (2004)

%   Author: Tanya Nair
%   Last Modified: March 23, 2016

%% Constants
% close all; clear all;
segType = {'otsu', 'kmeans'}; segIdx = 1;

%% Load brain
imDir = 'C:\Users\Lakshmi\OneDrive\Documents\ecse 626\project\tanyaCode\disc1_OAS1_0001_MR1_processed\MPRAGE\T88_111';
hdr = [imDir '\OAS1_0001_MR1_mpr_n4_anon_111_t88_masked_gfc.hdr'];
hdrInfo = analyze75info(hdr);
I = double(analyze75read( hdrInfo ));
% I_orig = double(X(:,:,100));
brainMask = I; brainMask(find(I==0))=1;

betascale = .055;
% BETA = 1.25*[1e-3 .2 .1];     % weighting of the V2 potentials for Beijing World Park 8.JPG
% BETA = betascale*[.5 .07 .07 .2];     % weighting of the V2 potentials

% BETA = betascale*[.75 1.13 1.2 .75];     % weighting of the V2 potentials for Beijing World Park 8.JPG
% BETA = betascale*[.8 1 1.1 .095];     % weighting of the V2 potentials for Beijing World Park 8.JPG
BETA = [0 .9 1.05 1];     % weights for neighbourhood potentials
ALPHA = [0 1.1 1 1];    % weights for unary potentials
NCOMPONENTS = 4;
MAXITER_EM = 100;
MAXITER_ICM = 10;
STOPPERCENT = .5;

% BETA = [1 1];
% ALPHA = 1.25*[1.5e-3 .2];

%% Load toy image
% MYIMG = 'Beijing World Park 8.JPG';
% NCOMPONENTS = 2;
% MAXITER_EM = 100;
% MAXITER_ICM = 10;
% STOPPERCENT = 5;
% 
% I = im2double(imread( MYIMG ));
% I_orig = im2double(rgb2gray(imread(MYIMG)));
% 
% betascale = .065;
% BETA = [.1 .1];     % weighting of the V2 potentials for Beijing World Park 8.JPG
% ALPHA = 1.25*[1 1];
%%
IMDIMS = size(I);
[I_initSeg, model] = getInitSeg( I, NCOMPONENTS, segType{segIdx}, brainMask );
labels = I_initSeg;

% [pxgn, labels, energy{t}] = runICM( I, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA );
%%
Umodel1 = cell(1,NCOMPONENTS);
for class=1:NCOMPONENTS
    Umodel1{class} = log(normpdf( I(:), model.mu(class), sqrt(model.sig(class))));
    Umodel1{class} = Umodel1{class}';
%     Umodel1{class} = normc(Umodel1{class}');
end
%%
Usum = zeros(1,MAXITER_ICM);
for iter=1:MAXITER_ICM
    % calculate energies
    U = cell(1, NCOMPONENTS);
    Umodel2 = calcEnergy2(labels, IMDIMS);
    Umodel2 = calcEnergy3D( labels, IMDIMS);
    for class=1:NCOMPONENTS
        U{class} = ALPHA(class)*Umodel1{class} + BETA(class)*(Umodel2);
    end
    [Umin, labels] = min(cell2mat(U),[],2);
%     labels(find(brainMask==1))=1;
    Usum(iter) = sum(Umin);
end
energy = sum(Umin);


%%
% IMDIMS = size(I_orig);
% I_orig = reshape(I_orig, [1 numel(I_orig)]);
% % the model contains the mu & sig
% [I_initSeg, model] = getInitSeg( I_orig, NCOMPONENTS, segType{segIdx}, brainMask);
% 
% [I_finalSeg, model, energy, emIter] = runHMRF( I_orig, I_initSeg, model, brainMask, NCOMPONENTS, ...
%                                 MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA, STOPPERCENT);
% 
% I_finalSeg(brainMask==1)=1;
%% display results
figure;
subplot(131); imagesc( reshape(I_orig, IMDIMS) ); title( 'original' );
subplot(132); imagesc( reshape(I_initSeg, IMDIMS) ); title( sprintf('initial seg: %s', segType{segIdx}) );
subplot(133); imagesc( reshape(I_finalSeg, IMDIMS) ); title( ...
                sprintf( 'final seg, %i EM, %i ICM, %.3f beta ', MAXITER_EM, MAXITER_ICM, betascale )); 
           
% figure;
% subplot(121); plot(energy); title( 'Energy vs. EM Iteration' ); hold on; plot( energy, 'ro');
% subplot(122); imagesc( reshape(I_finalSeg, IMDIMS) ); title( ...
%                 sprintf( 'final seg, %i EM, %i ICM, %.3f beta ', emIter, MAXITER_ICM, betascale ));

