function [Ty,diff,p]=yueneeg(a,b,percent,alpha)

% function [Ty,diff,p]=yueneeg(a,b,percent,alpha)
%
% Computes Ty (Yuen's T statistic).
% Ty=(tma-tmb) / sqrt(da+db), where tma & tmb are trimmed means of a & b,
% and da & db are yuen's estimate of the standard errors of tma & tmb.
% Ty is distributed approximately as Student's t with estimated degrees of freedom, df.
% The p-value (p) is the probability of obtaining a t-value whose absolute value
% is greater than Ty if the null hypothesis (i.e., tma-tmb = mu) is true.
% In other words, p is the p-value for a two-tailed test of H0: tma-tmb=0;
%
% a & b are matrices EEGLAB format with trials as the third dimension; 
% percent must be a number between 0 & 100.
%
%   Default values:
%   percent = 20;
%   alpha = 0.05; 
%
% See Wilcox (2005), Introduction to Robust Estimation and Hypothesis
% Testing (2nd Edition), page 159-161 for a description of the Yuen
% procedure.
% You can also check David C. Howell, Statistical Methods for Psychology,
% sixth edition, p.43, 323, 362-363.
%
% Original code by Prof. Patrick J. Bennett, McMaster University
% Added CI output, various editing, GAR, University of Glasgow, Dec 2007
% Vector version for EEG: GAR - University of Glasgow - Nov 2009
%
% See also YUEN

if nargin<4;alpha=.05;end
if nargin<3;percent=20;end

if isempty(a) || isempty(b) 
    error('yuen:InvalidInput', 'data vectors cannot have length=0');
end

if (percent >= 100) || (percent < 0)
    error('yuen:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

% sort
% boundaries, df
% winsorise, variance
% trim, mean
% ttest

% number of trials
na=size(a,3);
nb=size(b,3);

% number of items to winsorize and trim
ga=floor((percent/100)*na);
gb=floor((percent/100)*nb);

% winsorise a
asort=sort(a,3);
wa=asort;
wa(:,:,1:ga+1)=repmat(asort(:,:,ga+1),[1 1 ga+1]);
wa(:,:,na-ga:end)=repmat(asort(:,:,na-ga),[1 1 ga+1]);
wva=var(wa,0,3);

% winsorise b
bsort=sort(b,3);
wb=bsort;
wb(:,:,1:gb+1)=repmat(bsort(:,:,gb+1),[1 1 gb+1]);
wb(:,:,nb-gb:end)=repmat(bsort(:,:,nb-gb),[1 1 gb+1]);
wvb=var(wb,0,3);

% yuen's estimate of standard errors for a and b
ha=na-2*ga;
da=((na-1)*wva)/(ha*(ha-1));

hb=nb-2*gb;
db=((nb-1)*wvb)/(hb*(hb-1));

% trimmed means
ma=mean(asort(:,:,(ga+1):(na-ga)),3);
mb=mean(bsort(:,:,(gb+1):(nb-gb)),3);

diff=ma-mb;

Ty=diff ./ sqrt(da+db);

df= ((da+db).^2) ./ ( (da.^2/(ha-1)) + (db.^2/(hb-1)) );

p=2*(1-tcdf(abs(Ty),df)); % 2-tailed probability

% t=tinv(1-alpha./2,df); % 1-alpha/2 quantile of Student's distribution with df degrees of freedom
% CI(1)=(ma-mb)-t.*sqrt(da+db); 
% CI(2)=(ma-mb)+t.*sqrt(da+db);

return
