function out = charpoly4x4_sym(z)

% compute coefficients of characteristic polynomial of hermitian 4x4
% matrices

% !!!! no check of symmetry is performed


out = zeros(size(z,1),5);

z_abssq = abs(z).^2;

tr1 = trace4x4_1(z);
tr2 = trace4x4_2(z_abssq);
tr3 = trace4x4_3(z, z_abssq);

out(:,1) = 1;
out(:,2) = -tr1; 
out(:,3) = 0.5.*(tr1.^2-tr2);
out(:,4) = -(tr1.^3)./6 + (tr1.*tr2)./2 - tr3./3;
out(:,5) = real(det4x4(permute(z,[2 3 1])));

function tr1 = trace4x4_1(z)

tr1 = z(:,1,1)+z(:,2,2)+z(:,3,3)+z(:,4,4);

function tr2 = trace4x4_2(z_abssq)

tr2 = sum(sum(z_abssq,3),2);

function tr3 = trace4x4_3(z, z_abssq)

tr3 =  (z(:,1,1).^3+z(:,2,2).^3+z(:,3,3).^3+z(:,4,4).^3) + ...
      3.*(z(:,1,1).*sum(z_abssq(:,[2 3 4],1),2) + z(:,2,2).*sum(z_abssq(:,[1 3 4],2),2) + ...
          z(:,3,3).*sum(z_abssq(:,[1 2 4],3),2) + z(:,4,4).*sum(z_abssq(:,[1 2 3],4),2)) + ...
          6.*real(z(:,3,1).*z(:,1,2).*z(:,2,3) + z(:,4,1).*z(:,1,2).*z(:,2,4) + ...
                  z(:,4,1).*z(:,1,3).*z(:,3,4) + z(:,4,2).*z(:,2,3).*z(:,3,4));



