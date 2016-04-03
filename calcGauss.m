function gauss = calcGauss( x, mean, var)

diff = abs(x - mean);
gauss = diff.*diff - (2*var) + log(sqrt(2*pi*var)) ;
gauss = gauss';
gauss = gauss/max(gauss,[],1);


gauss2 = 1/(sqrt(var*2*pi)) * exp(- (x-mean).*(x-mean) / (2*var));