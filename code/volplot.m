function [dat, map] = volplot(dat, sel, cscale, cmap)

% VOLPLOT make 2D or 3D plot of volumetric data (e.g. MRI)
% that is defined on a regular orthogonal grid
% 
% volplot(dat, sel) or
% volplot(x, y, z, dat, sel)
% volplot(x, y, z, dat, sel, caxis)
%
% where sel is one of
%   [x, y, z]     intersection through the three orthogonal directions
%   index         linear index of the voxel of interest
%   'min'         intersection at the minimum
%   'max'         intersection at the maximum
%   'center'      intersect at the center of each axis
%   'interactive' intersect at the center, then go into interactive mode
%   'maxproject'  project the maximum value along each orthogonal direction
%   'sumproject'  integrated value along each orthogonal direction (glassbrain)
%   'montage'     show all slices
% and caxis is the [min max] used for the color scaling
% 
% See also TRIPLOT, LINEPLOT (in ~roberto/matlab/misc)
% See also NDGRID

% Copyright (C) 2003, Robert Oostenveld
% 
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: volplot.m 6314 2012-08-03 12:24:16Z jansch $

if nargin<2 || isempty(sel)
  sel = 'interactive'; 
end
if nargin<3 || isempty(cscale)
  cscale = [];
end
if nargin<4 || isempty(cmap)
  cmap = 'jet';
end

x = 1:size(dat,1);
y = 1:size(dat,2);
z = 1:size(dat,3);

if ~isempty(cscale)
  cmin = cscale(1);
  cmax = cscale(2);
else
  % ensure same color scaling for all figures
  cmin = min(dat(:));
  cmax = max(dat(:));
end
  
% convert from 1-d to 3-d array
if any(size(dat)==1)
  dim(1) = length(x);
  dim(2) = length(y);
  dim(3) = length(z);
  dat = reshape(dat, dim);
else
  dim = size(dat);
end

% convert the selection to the indices of the x/y/z intersection 
if ischar(sel) && strcmp(sel, 'min')
  [minval, minindx] = min(dat(:));
  [xi, yi, zi] = ind2sub(size(dat(:,:,:,1,1)), minindx);
elseif ischar(sel) && strcmp(sel, 'max')
  [maxval, maxindx] = max(dat(:));
  [xi, yi, zi] = ind2sub(size(dat(:,:,:,1,1)), maxindx);
elseif ischar(sel) && strcmp(sel, 'center')
  xi = round(length(x)/2);
  yi = round(length(y)/2);
  zi = round(length(z)/2);
elseif ischar(sel) && strcmp(sel, 'interactive')
  xi = round(length(x)/2);
  yi = round(length(y)/2);
  zi = round(length(z)/2);
  fc1 = 1;
  fc2 = 1;
elseif ~ischar(sel) && length(sel)==1
  [xi, yi, zi] = ind2sub(dim(1:3), sel);
  fc1 = 1;
  fc2 = 1;
  
else
  xi = nearest(x, sel(1));
  yi = nearest(y, sel(2));
  zi = nearest(z, sel(3));
  if numel(sel)>=5
    fc1 = sel(4);
    fc2 = sel(5);
  else
    fc1 = 1;
    fc2 = 1;
  end
end

cscaleorig = cscale;

% start the plotting
if strcmp(sel, 'interactive')
  xc = x(xi);
  yc = y(yi);
  zc = z(zi);
  
  set(gcf, 'pointer', 'cross');
  
  nas = [];
  lpa = [];
  rpa = [];
  while(1)
    fprintf('============================================================\n');
    fprintf('click with mouse button to reslice the display to a new position\n');
    fprintf('press n/l/r on keyboard to record a fiducial position\n');
    fprintf('press q on keyboard to quit interactive mode\n');
    volplot(dat, [xc yc zc fc1 fc2], cscale, cmap);
    drawnow;
    
    try
      waitforbuttonpress;
      ax = gca;
      p  = get(ax, 'currentpoint');
      
      d1 = p(1,1);
      d2 = p(1,2);
      
      key = get(gcf, 'currentcharacter');
      
      %fprintf('d1=%d, d2=%d, ax=%s,%s\n', d1, d2, get(get(ax,'xlabel'),'string'), get(get(ax,'ylabel'),'string'));
    catch 
      key='q'; 
    end
    
    if key=='q'
      % rename the first output argument of this function
      h = nas;
      break;
    elseif key=='l'
      lpa = [xc yc zc];
    elseif key=='r'
      rpa = [xc yc zc];
    elseif key=='n'
      nas = [xc yc zc];
    elseif key=='x'
      selx= round((xc-0):(xc+0));
      sely= round((yc-0):(yc+0));
      selz= round((zc-0):(zc+0));
      dat(selx,sely,selz) = 0;
    elseif key=='X'
      selx= round((xc-1):(xc+1));
      sely= round((yc-1):(yc+1));
      selz= round((zc-1):(zc+1));
      dat(selx,sely,selz) = 0;
    elseif key==3
      % right mouse
      switch l1
        case {'x' 'y' 'z'}
          cscale(1) = min(reshape(dat(:,:,:,fc1,fc2),[],1));
          cscale(2) = max(reshape(dat(:,:,:,fc1,fc2),[],1));
        case {'freq1' 'freq2'}
          cscale(1) = min(reshape(dat(xi,yi,zi,:,:),[],1));
          cscale(2) = max(reshape(dat(xi,yi,zi,:,:),[],1));
      end
      
    else         
      % update the view to a new position
      l1 = get(get(ax, 'xlabel'), 'string');
      l2 = get(get(ax, 'ylabel'), 'string');
      switch l1
        case 'x'
          xc = d1;
          %cscale = cscaleorig;
        case 'y'
          yc = d1;
          %cscale = cscaleorig;
        case 'z'
          zc = d1;
          %cscale = cscaleorig;
        case 'freq1'
          fc1 = round(d2);
          fc2 = round(d1);
          
        case 'freq2'
          fc1 = round(d2);
          fc2 = round(d1);
          %cscale(1) = min(reshape(dat(:,:,:,fc1,fc2),[],1));
          %cscale(2) = max(reshape(dat(:,:,:,fc1,fc2),[],1));
      end
      switch l2
        case 'x'
          xc = d2;
          %cscale = cscaleorig;
        case 'y'
          yc = d2;
          %cscale = cscaleorig;
        case 'z'
          zc = d2;
          %cscale = cscaleorig;
      end
    end
    if ~isempty(nas), fprintf('nas = [%f %f %f]\n', nas); else fprintf('nas = undefined\n'); end 
    if ~isempty(lpa), fprintf('lpa = [%f %f %f]\n', lpa); else fprintf('lpa = undefined\n'); end 
    if ~isempty(rpa), fprintf('rpa = [%f %f %f]\n', rpa); else fprintf('rpa = undefined\n'); end 
  end
  
elseif strcmp(sel, 'montage')
  % make plot of x-y slices for all z values
  maxval = max(dat(:));
%  for z=1:size(dat,3)
%    % convert to 4D image for montage display
%    % transpose to correct for x-y axis change in Matlab image function
%    img(:,:,1,z) = transpose(dat(:,:,z));
%  end
%  montage(img);
%  axis xy

  % this works if java fails
  siz    = size(dat);
  ndiv   = [ceil(sqrt(siz(3))) floor(sqrt(siz(3)))];
  map    = ones(siz(2)*ndiv(2),siz(1)*ndiv(1))*dat(1);
  for k=1:siz(3)
    [nx,ny] = ind2sub(ndiv,k);
    map(siz(2)*(ny-1)+1:siz(2)*ny,siz(1)*(nx-1)+1:siz(1)*nx) = dat(:,:,k)';
  end
  h = imagesc(map);
  set(h, 'alphadata', ~isnan(map))  % make NaN values behave
  axis xy;axis equal;axis off;
  caxis([cmin cmax]);
  colormap(cmap)

elseif strcmp(sel, 'sumproject')
  % make plot of integrated-value projection along the thee orthogonal directions
  delete(subplot(2,2,4));   % delete the old colorbar
  h1 = subplot(2,2,1);
  h2 = subplot(2,2,2);
  h3 = subplot(2,2,3);
  h4 = subplot(2,2,4);  % this will  be the new colorbar

  % change not-a-number values to zero
  dat(find(isnan(dat(:)))) = 0;

  siz  = size(dat);
  
  if isempty(cscale)
    cmin = 0;
    cmax = max(max(max(reshape(sum(dat,1),siz(2)*siz(3),1)),max(reshape(sum(dat,2),siz(1)*siz(3),1))),max(reshape(sum(dat,3),siz(1)*siz(2),1)));
  else
    cmin = cscale(1);
    cmax = cscale(2);
  end
  
  subplot(h1);
  imagesc(x, z, squeeze(sum(dat, 2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h2);
  imagesc(y, z, squeeze(sum(dat, 1))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('y'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h3);
  imagesc(x, y, squeeze(sum(dat, 3))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('y');
  caxis([cmin cmax]);

  subplot(h4);
  colorbar;
  caxis([cmin cmax]);
  xlabel('colorscale')

elseif strcmp(sel, 'maxproject')
  % make plot of maximum-value projection along the thee orthogonal directions
  delete(subplot(2,2,4));   % delete the old colorbar
  h1 = subplot(2,2,1);
  h2 = subplot(2,2,2);
  h3 = subplot(2,2,3);
  h4 = subplot(2,2,4);  % this will  be the new colorbar

  subplot(h1);
  imagesc(x, z, squeeze(max(dat, [], 2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h2);
  imagesc(y, z, squeeze(max(dat, [], 1))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('y'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h3);
  imagesc(x, y, squeeze(max(dat, [], 3))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('y');
  caxis([cmin cmax]);

  subplot(h4);
  colorbar(h4, 'peer', h1);
  xlabel('colorscale')

else
  % make plot of three orthogonal slices intersecting at [xi yi zi]
  if ~exist('xi', 'var') || ~exist('yi', 'var') || ~exist('zi', 'var')
    error('nothing to plot, no selection given')
  end

  fprintf('value of %f in voxel %d at [%.02f %.02f %.02f]\n', double(dat(xi, yi, zi)), sub2ind(dim, xi, yi, zi), x(xi), y(yi), z(zi));

  warning off
  delete(subplot(2,2,4));   % delete the old colorbar
  warning on
  h1 = subplot(2,2,1);
  h2 = subplot(2,2,2);
  h3 = subplot(2,2,3);
  h4 = subplot(2,2,4);  % this will  be the new colorbar

  subplot(h1);
  imagesc(x, z, squeeze(dat(:,yi,:,fc1,fc2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('z');
  caxis([cmin cmax]);
  crosshair([x(xi) z(zi)], 'color', 'yellow');
  colormap(cmap);
  
  subplot(h2);
  imagesc(y, z, squeeze(dat(xi,:,:,fc1,fc2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('y'); ylabel('z');
  caxis([cmin cmax]);
  crosshair([y(yi) z(zi)], 'color', 'yellow');
  colormap(cmap);
  
  subplot(h3);
  imagesc(x, y, squeeze(dat(:,:,zi,fc1,fc2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('y');
  caxis([cmin cmax]);
  crosshair([x(xi) y(yi)], 'color', 'yellow');
  colormap(cmap);
  
  if ndims(dat)>3
    subplot(h4);
    imagesc(squeeze(dat(xi,yi,zi,:,:)));
    xlabel('freq1'); ylabel('freq2');
    caxis([cmin cmax]);
  end
end

