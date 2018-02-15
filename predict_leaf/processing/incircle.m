function [C,R] = incircle(x,y)
if (nargin < 1) || (nargin > 2)
  error('INCIRCLE:improperarguments', ...
    'incircle requires exactly 1 or 2 arguments')
elseif (nargin == 2) && isvector(x) && isvector(y) && (numel(x) == numel(y))
  % a pair of vectors
  xy = [x(:),y(:)];
  % compute the hull. I prefer convhulln.
  edges = convhulln(xy);
elseif (nargin == 1) && (size(x,2) == 2)
  % a single list of points as rows of x
  xy = x;
  edges = convhulln(xy);
elseif (nargin == 2) && (size(x,2) == 2) && (size(y,2) == 2)
  % y must be a list of edges from convhulln
  xy = x;
  edges = y;
elseif (nargin == 2) && (size(x,2) == 2) && isvector(y) && (y(1) == y(end))
  % y must be a list of edges from convhull
  xy = x;
  I = y(:);
  % in case the first point was not wrapped in the polygon
  if I(1) ~= I(end)
    I = [I;I(1)];
  end
  edges = [I(1:(end-1)),I(2:end)];
else
  % none of the forms I allow seem to fit.
  % give up and throw an error.
  error('INCIRCLE:invaliddata', ...
    'x and y do not seem to fit into any of the allowed forms for incircle input')
end
ne = size(edges,1);

% the end points of each edge are...
A = xy(edges(:,1),:);
B = xy(edges(:,2),:);

% the normal vector to each edge
N = (B - A)*[0 1;-1 0];

% normalize to unit length
L = sqrt(sum(N.^2,2));
N = N./[L,L];

% a central point inside the hull itself
C0 = mean(A,1);

k = sum(N.*bsxfun(@minus,C0,A),2) < 0;
N(k,:) = -N(k,:);

Aeq = [N,zeros(ne,1),-eye(ne)];
beq = sum(N.*A,2);

% lower bounds only for the slack variables
LB = [-inf; -inf; 0; zeros(ne,1)];
% there are no upper bounds
UB = [];

% inequalities defined by the slack variable
% constraints
A = [zeros(ne,2),ones(ne,1),-eye(ne)];
b = zeros(ne,1);

% the objective is just -R
f = zeros(ne+3,1);
f(3) = -1;

% just turn off the output message
options = optimset('linprog');
options.Display = 'off';

% linprog does all of the hard work.
result = linprog(f,A,b,Aeq,beq,LB,UB,[],options);

% unpack the circle parameters
C = result(1:2)';
R = result(3);



