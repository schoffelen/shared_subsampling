load leadfield
L  = cat(2,leadfield.leadfield{:});
L1 = L(:,1:2:end);
L2 = L(:,2:2:end);

cd ~/matlab/misc
C = vol2connmat(leadfield,6);
C = C(leadfield.inside,leadfield.inside);

inside = find(leadfield.inside);

% this is good enough to get the maximum correlation of neighbour leadfield
% components (although the optimally aligned orientation is not taken into
% account, who cares?
for k = 1:4416
  c = corr(L1(:,[k find(C(k,:))]));  
  cmax1(k,1) = max(abs(c(2:end,1)));
  c = corr(L2(:,[k find(C(k,:))]));
  cmax2(k,1) = max(abs(c(2:end,1)));
  c = corr([L1(:,k) L2(:,[find(C(k,:))])]);
  cmax3(k,1) = max(abs(c(2:end,1)));
  c = corr([L2(:,k) L1(:,[find(C(k,:))])]);
  cmax4(k,1) = max(abs(c(2:end,1)));
end
cmax = max(cmax1,cmax2,cmax3,cmax4);
figure;histogram(cmax);
xlabel('maximal correlation leadfield component');

% alternatively, get the maximum condition number of the 4x4 matrices
for k = 1:4416
  sel = find(C(k,:));
  condmax(k,1) = 0;
  for m = 1:numel(sel)
    tmplf = [L1(:,[k sel(m)]) L2(:,[k sel(m)])];
    condmax(k) = max(condmax(k), cond(tmplf'*tmplf));
  end
end
figure;histogram(log10(dum(inside))); xlabel('log_1_0 condition number'); ylabel('# dipoles');

