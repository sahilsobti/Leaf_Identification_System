function [J, grad] = lrCostFunction(theta, X, y, lambda)

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));
h=sigmoid(X * theta);
ntheta=theta(2:end,1);
J=(1/m).*((-y'*log(h))-((1-y)'*log(1-h)))+(lambda/(2*m)).*(ntheta' * ntheta);
grad1 = (1/m) .* ((h-y)' * X(:,1));
grad2 = ((1/m) .* ((h-y)'* X(:,2:end))) + ((lambda/m) .* theta(2:end,1)');
grad = [grad1 grad2];
% =============================================================
grad = grad(:);

end
