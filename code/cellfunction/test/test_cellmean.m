function test_cellmean

% x as a 1xn cell-array
x = {rand(3,4) rand(3,4)};

m1 = mean(cat(1, x{:}));
m2 = mean(cat(2, x{:}),2);

m1c = cellmean(x, 1);
m2c = cellmean(x, 2);

assert(norm(m1-m1c)./norm(m1)<eps);
assert(norm(m2-m2c)./norm(m2)<eps);

% x as a nx1 cell-array
x = x(:);

m1c = cellmean(x, 1);
m2c = cellmean(x, 2);

assert(norm(m1-m1c)./norm(m1)<eps);
assert(norm(m2-m2c)./norm(m2)<eps);
