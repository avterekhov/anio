function f = polyfunc(x, w)
% f = polyfunc(x, w)
%
% Additive polynomial function defined by its coefficients w.
%
% f = w(:,1).*x + w(:,2).*x.^2 + ... + w(:,n).*x.^n
%
% Copyright (C) 2011, Alexander V. Terekhov, avterekhov@gmail.com
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301 USA.

pmax = size(w,2)-1;
f = zeros(size(x));
for p = 0:pmax
  f = f + repmat(w(:,p+1),1,size(x,2)).*x.^p;
end;
