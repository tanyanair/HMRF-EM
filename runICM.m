function [pxgn, labels, energy] = runICM( I_orig, labels, model, brainMask, NCOMPONENTS, MAXITER_ICM, IMDIMS, BETA, ALPHA )

brainMask = reshape(brainMask, [1 numel(brainMask)]);

% calculate the likelihood for the data with the updated model parameters
Umodel1 = cell(1,NCOMPONENTS);
for class=1:NCOMPONENTS
%     Umodel1{class} = calcGauss( I_orig, model.mu(class), model.sig(class));
    Umodel1{class} = log(normpdf( I_orig(:), model.mu(class), sqrt(model.sig(class))));
    Umodel1{class} = Umodel1{class}';
%     Umodel1{class} = normc(Umodel1{class}');
end

Usum = zeros(1,MAXITER_ICM);
for iter=1:MAXITER_ICM
%    fprintf( '     ICM Iteration: %i\n', iter);
    % calculate energies
    U = cell(1, NCOMPONENTS);
    Umodel2 = calcEnergy2(labels, IMDIMS, NCOMPONENTS);
    for class=1:NCOMPONENTS
        U{class} = ALPHA(class)*Umodel1{class} + BETA(class)*(Umodel2);
    end
    [Umin, labels] = min(cell2mat(U),[],2);
    labels(brainMask==1)=1;
    Usum(iter) = sum(Umin);
end
energy = sum(Umin);

% P( x_i=label | neighbourhood of x_i )
pxgn = cell(1,NCOMPONENTS);
for class = 1:NCOMPONENTS
   pxgn{class} = normc((exp(U{class})));
end

% figure;
% plot(Usum); 
% drawnow;
