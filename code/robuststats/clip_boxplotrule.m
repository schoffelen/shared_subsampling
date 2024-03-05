function out = clip_boxplotrule(in)

[m,n]    = size(in);

[t1, t2] = percthreshold(in, 0.25, 2);
iq       = t2-t1;
upper    = t2+1.5*iq;
%lower    = t1-1.5*iq; % don't clip the lower

out      = min(in, upper*ones(1,n));

