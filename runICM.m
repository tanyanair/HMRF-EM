%% runICM obtains MRF-MAP labels for the image by minimizing the posterior 
%   energy of the image labels.
%
% Author: Tanya Nair
% Last Modified: May 7, 2016

function [pxgn, labels, energy] = runICM( I, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA )

% calculate the likelihood for the data with the updated model parameters
% (this is prior)
Umodel1 = cell(1,NCOMPONENTS);
for class=1:NCOMPONENTS
    Umodel1{class} = -log(normpdf( I(:), model.mu(class), sqrt(model.sig(class))));
end

Usum = zeros(1,MAXITER_ICM);
for iter=1:MAXITER_ICM
%    fprintf( '     ICM Iteration: %i\n', iter);
    % calculate energies
    U = cell(1, NCOMPONENTS);
    Umodel2 = calcEnergy3D(labels, brainMask, IMDIMS, NCOMPONENTS, BETA);
    for class=1:NCOMPONENTS
        U{class} = ALPHA(class)*Umodel1{class} + Umodel2{class};
    end
    [Umin, labels] = min(cell2mat(U),[],2);
    labels(brainMask==1)=0;
    Usum(iter) = sum(Umin);
end
energy = sum(Umin);

% P( x_i=label | neighbourhood of x_i )
pxgn = cell(1,NCOMPONENTS);
for class = 1:NCOMPONENTS
   pxgn{class} = normr((exp(-U{class})));
end

% figure;
% plot(Usum); title( 'Energy Sum per ICM Iteration' );
% xlabel( 'ICM Iteration' ); ylabel( 'Sum of U' );
% drawnow;
