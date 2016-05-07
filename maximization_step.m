% maximization_step performs maximum likelihood calculations for each 
%   class. This is the "M" step in EM.
%
% Author: Tanya Nair
% Last Modified: May 7, 2016

function [model, logli_img] = maximization_step( img, pxgn, model, NCOMPONENTS )

img = img(:);
% calculate plgy = P( label | image );
gauss = cell( 1,NCOMPONENTS );
plgy = cell( 1, NCOMPONENTS );
p_y = 0;
for k=1:NCOMPONENTS
   gauss{k} = (normpdf( img, model.mu(k), sqrt(model.sig(k)) ));
   p_y = p_y + gauss{k};
end
for k=1:NCOMPONENTS
  plgy{k} = ((pxgn{k}.*gauss{k}/sum(p_y)));
end

% update model parameters
for k=1:NCOMPONENTS
    model.mu(k) = sum(img.*plgy{k})/sum(plgy{k});
    model.sig(k) = sum((plgy{k}).*(img - model.mu(k)).^2)/sum(plgy{k});
end

logli_img = p_y;