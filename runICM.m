function [pxgn, labels, energy] = runICM( I, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA )

% calculate the likelihood for the data with the updated model parameters
Umodel1 = cell(1,NCOMPONENTS);
for class=1:NCOMPONENTS
    Umodel1{class} = log(normpdf( I(:), model.mu(class), sqrt(model.sig(class))));
end

Usum = zeros(1,MAXITER_ICM);
for iter=1:MAXITER_ICM
%    fprintf( '     ICM Iteration: %i\n', iter);
    % calculate energies
    U = cell(1, NCOMPONENTS);
    Umodel2 = calcEnergy3D(labels, brainMask, IMDIMS, NCOMPONENTS);
    for class=1:NCOMPONENTS
        U{class} = ALPHA(class)*Umodel1{class} + BETA(class)*(Umodel2);
    end
    [Umin, labels] = min(cell2mat(U),[],2);
    labels(brainMask==1)=0;
    Usum(iter) = sum(Umin);
end
energy = sum(Umin);

% P( x_i=label | neighbourhood of x_i )
pxgn = cell(1,NCOMPONENTS);
for class = 1:NCOMPONENTS
%    pxgn{class} = normc((exp(U{class})));
    pxgn{class} = normr(U{class});
end

% figure;
% plot(Usum); 
% drawnow;
