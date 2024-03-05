function d = inv4x4(x)

% INV4X4 computes inverse of matrix x, using explicit analytic definition
% if size(x) = [4 4 K M]

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
% $Id$

siz = size(x);
if all(siz(1:2)==4)
  adjx = [ det3x3(x([2 3 4],[2 3 4],:,:)) -det3x3(x([1 3 4],[2 3 4],:,:))  det3x3(x([1 2 4],[2 3 4],:,:))  -det3x3(x([1 2 3],[2 3 4],:,:)); ...
          -det3x3(x([2 3 4],[1 3 4],:,:))  det3x3(x([1 3 4],[1 3 4],:,:)) -det3x3(x([1 2 4],[1 3 4],:,:))   det3x3(x([1 2 3],[1 3 4],:,:)); ...
	         det3x3(x([2 3 4],[1 2 4],:,:)) -det3x3(x([1 3 4],[1 2 4],:,:))  det3x3(x([1 2 4],[1 2 4],:,:))  -det3x3(x([1 2 3],[1 2 4],:,:)); ...
          -det3x3(x([2 3 4],[1 2 3],:,:))  det3x3(x([1 3 4],[1 2 3],:,:)) -det3x3(x([1 2 4],[1 2 3],:,:))   det3x3(x([1 2 3],[1 2 3],:,:)); ...
           ];
  denom = det4x4(x);
  d     = adjx./denom([1 1 1 1],[1 1 1 1],:,:);
else
  error('cannot compute slicewise inverse');
  % write for loop for the higher dimensions, using normal inv
end
