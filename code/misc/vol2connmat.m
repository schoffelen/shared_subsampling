function [connmat] = vol2connmat(source, conndef)

% VOL2CONNMAT creates the matrix containing the voxel
% neighbourhood structure.
% Use as [connmat] = vol2connmat(source, conndef)
% source should contain the fields dim and inside
% supported values for conndef are 6 (default), 18 and 26
% otherwise conndef is interpreted as a scalar distance measure.
% in that case a matrix is returned containing the distance of 
% those voxels to the reference voxel. only voxels with a distance 
% (in source-units) <= conndef to the reference voxel is taken into account.

flag = 0;
if nargin==1,                     conndef = 6; end
if ~ismember(conndef, [6 18 26]), flag    = 1; end

%---create kernel
if conndef==6
  kernel = false(3,3,3);
  kernel(2,2,1) = true;
  kernel(2,2,3) = true;
  kernel(2,1,2) = true;
  kernel(2,3,2) = true;
  kernel(1,2,2) = true;
  kernel(3,2,2) = true;
elseif conndef==18
  kernel = true(3,3,3);
  kernel(2,2,2) = false;
  kernel(3,3,3) = false;
  kernel(1,3,3) = false;
  kernel(3,1,3) = false;
  kernel(3,3,1) = false;
  kernel(1,1,3) = false;
  kernel(1,3,1) = false;
  kernel(3,1,1) = false;
  kernel(1,1,1) = false;
elseif conndef==26
  kernel = true(3,3,3);
  kernel(2,2,2) = false;
else
  dpos    = sqrt(sum((source.pos(2,:)-source.pos(1,:)).^2));
  nvox    = 2*ceil(conndef./dpos) + 1;
  kernel  = zeros(nvox, nvox, nvox);
  ind     = [1:nvox] - (nvox+1)/2;
  [x,y,z] = ndgrid(ind, ind, ind);
  kernel(:)  = dpos.*sqrt(x(:).^2 + y(:).^2 + z(:).^2);
  kernel(kernel>conndef) = 0;
end

%---allocate memory
if islogical(source.inside)
  source.inside = find(source.inside);
end
nvox    = prod(source.dim);
ninside = length(source.inside);

if ~flag
  connmat = logical(spalloc(nvox, nvox, sum(kernel(:)).*ninside));
  %---fill matrix
  mask = false(source.dim);
  mask(source.inside) = true;
  for j = source.inside(:)'
    tmpvol        = false(source.dim+2);
    %tmpvol(j)     = true;
    %connmat(:, j) = reshape(mask & imdilate(tmpvol, kernel), [nvox 1]);
    [x,y,z]       = ind2sub(source.dim, j);
    tmpvol(x:x+2,y:y+2,z:z+2) = kernel;
    connmat(:, j) = reshape(mask & tmpvol(2:end-1,2:end-1,2:end-1), [nvox 1]);
  end
else
  connmat = spalloc(nvox, nvox, sum(find(kernel(:))).*ninside);
  %---fill matrix
  mask = zeros(source.dim);
  mask(source.inside) = 1;
  siz    = size(kernel);
  offset = (siz(1)-1)/2;
  for j = source.inside(:)'
    tmpvol        = zeros(source.dim+2*offset);
    [x,y,z]       = ind2sub(source.dim, j);
    tmpvol(x:x+siz(1)-1, y:y+siz(2)-1, z:z+siz(3)-1) = kernel;
    connmat(:, j) = reshape(mask.*tmpvol(offset+1:end-offset,offset+1:end-offset,offset+1:end-offset), [nvox 1]);
  end
end
