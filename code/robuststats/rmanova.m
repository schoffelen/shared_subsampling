function [F, df1, df2] = rmanova(x, percent)

% RMANOVA(X, PERCENT) performs a repeated-measures anova Omnibus test
% on the trimmed means, according to Wilcox section 8.1
%
% Input arguments: 
%   x = data matrix N observations x M repetitions
%   percent = the amount of trimming

if nargin<2, 
  percent = 0.2;
end

% prepare some stuff
n = size(x, 1);
J = size(x, 2);
g = floor(n*percent);
h = n-2*g;

% get the mean of the original groups + the grand mean
mu_x    = trimmean(x, percent, 1);
mu_xall = mean(mu_x);

% winsorize the sample
y = winsorize(x, percent, 1);

% get some quantities from the winsorized groups
y_dot_j = mean(y, 1);
y_i_dot = mean(y, 2);
y_dot   = mean(y(:));

% get some sums-of-squares
Qc = h.*sum((mu_x - mu_xall).^2);
Qe = sum(sum( (y - ones(n,1)*y_dot_j - y_i_dot*ones(1,J) + y_dot).^2) );

F = ( Qc./(J-1) )./( Qe./((h-1)*(J-1)) ); 

% get the degrees of freedom
v       = cov(y);
v_dot   = sum(v(:))./(J^2);
v_d     = trace(v)./J;
v_j_dot = mean(v, 2);

A =  ((J^2)/(J-1))*((v_d-v_dot)^2);
B =  sum(v(:).^2) - 2*J*sum(v_j_dot.^2) + (J*v_dot)^2;

eps_hat   = A./B;
eps_tilde = (n*(J-1)*eps_hat-2) / ((J-1)*(n-1-(J-1)*eps_hat)); 

df1 = (J-1)*eps_tilde;
df2 = (J-1)*(h-1)*eps_tilde;
