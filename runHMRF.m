%% runHMRF runs the EM routine described in Zhang et al., 2001
%   The function iteratively runs the ICM algorithm to obtain MRF-MAP
%   labels, and then the maximization_step to obtain maximum likelihood
%   parameter estimates for the model.
%
% Author: Tanya Nair
% Last Modified: May 7, 2016

function [final_seg, model, energy, t] = runHMRF( I_orig, labels, model, brainMask, NCOMPONENTS, ...
                                MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA)                         
                            
energy = cell(1,MAXITER_EM);
final_seg = reshape(labels, IMDIMS);
for t=2:MAXITER_EM
    % expectation step: best guess of tissue labels (MAP-MRF with ICM)
    [pxgn, labels, energy{t}] = runICM( I_orig, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA );
    
    % do maximum likelihood estimation for the new model parameters
    [model, ~] = maximization_step( I_orig, pxgn, model, NCOMPONENTS);
   
    fprintf( 'EM Iter %i complete, Energy: %.1d \n', t-1, energy{t} );
    
    % check for convergence..
    if( energy{t} > energy{t-1} )
        break;
    end
    if( energy{t} < 0 )
        break;
    end
    final_seg = reshape( labels, IMDIMS );
    
end

% plot the energy over all the ICM iterations across all the EM iterations
energy = cell2mat(energy);
