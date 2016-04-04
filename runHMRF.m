function [final_seg, model, energy, t] = runHMRF( I_orig, labels, model, brainMask, NCOMPONENTS, ...
                                MAXITER_EM, MAXITER_ICM, IMDIMS, BETA, ALPHA, STOPPERCENT)

% stopCriteriaPercent = S;
percent_diff = 100;
logli_img = 100;
energy = cell(1,MAXITER_EM);
for t=1:MAXITER_EM
    % check for convergence..
    if( abs(percent_diff) < STOPPERCENT )
       break;
    end
    logli_prev = logli_img;

    % expectation step: best guess of tissue labels (MAP-MRF with ICM)
    [pxgn, labels, energy{t}] = runICM( I_orig, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA );
    
    % do maximum likelihood estimation for the new model parameters
    [model, logli_img] = maximization_step( I_orig, pxgn, model, NCOMPONENTS);
   
    percent_diff = -100*mean(( (logli_img - logli_prev)./(logli_prev) ));
%     fprintf( 'EM Iter %i complete, Logli %.1f, Change Logli %.2f%%\n', t, sum(logli_img), percent_diff);
    fprintf( 'EM Iter %i complete, Change Logli %.2f%%\n', t, percent_diff);
end
final_seg = reshape(labels,IMDIMS);

% plot the energy over all the ICM iterations across all the EM iterations
energy = cell2mat(energy);
% figure;
% plot(energy, 'r'); hold on; plot( energy, 'bo' );
% title( 'Energy vs. ICM Iteration');
