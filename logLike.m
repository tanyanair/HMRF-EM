function pixelScores = logLike( model, X, ncomponents )

% from slide 29 of the GMM-EM slides
%  p (x | theta ) = product_i_pixels[ sum_k_components [ a_k * N(x_i; mu_i, sigma_i) ] ]
pixelScores = zeros(size(X,1), size(X,2));
c=.01;
for k=1:ncomponents
    try
        pixelScores = log(model.alpha{k}*normpdf(X, model.mu{k}, sqrt(model.sig{k})));
    catch
        warning( 'adding regularization term to sigma because its too small');
        pixelScores = log(model.alpha{k}*mvnpdf(X, model.mu{k}, abs(model.sig{k}+c)));
    end
end