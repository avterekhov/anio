function w = anio(X, C, p, ffunc)
% w = anio(X, C)
%
% Approximates the cost function with the 2nd order additive polynomial from
% the data X obtained under constraints
%
%   C x = b
%
% The cost function has the form:
%
% J(x(1), ..., x(n)) = w(1,1) x(1)   + ... + w(n,1) x(n) +
%             + 1/2 * (w(1,2) x(1)^2 + ... + w(n,2) x(n)^2)
%
% w = anio(X, C, p)
%
% Approximates the cost function with the p-th order
% polynomials. The function has the form
%
% J(x(1), ..., x(n)) = w(1,1) x(1)   + ... + w(n,1) x(n)    +
%             + 1/2 * (w(1,2) x(1)^2 + ... + w(n,2) x(n)^2) +
%             .............................................
%             + 1/p * (w(1,p) x(1)^p + ... + w(n,p) x(n)^p)
%
%
% w = anio(X, C, p, func)
% Approximates the cost function with by an additive function,
% whose gradient is defined by func(x,w).
%
% The approximation method is described in:
%
% Terekhov, A. V., & Zatsiorsky, V. M. (2011). Analytical and
% numerical analysis of inverse optimization problems: conditions of
% uniqueness and computational methods. Biological cybernetics,
% 104(1-2), 75-93.
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


if nargin < 4
  ffunc = @(x,w)polyfunc(x,w);
end;

if nargin < 3
  p = 2;
end;


n = size(C,2);
k = size(C,1);
N = size(X,1);

Cbrev = eye(size(C,2))-C'*inv(C*C')*C;

A = [zeros(1,n), ones(1,n*(p-1));
     eye(n)-Cbrev, zeros(n,n*(p-1))];
b = [1; zeros(n,1)];

A = [zeros(1,n), ones(1,n*(p-1));
     eye(n)-Cbrev, zeros(n,n*(p-1))];
A = fullrank(A);
b = [1; zeros(k,1)];

X2 = zeros(n*p,n*p);
for s = 1:N
  Xs = zeros(n,n*p);
  Hs = zeros(n,p);
  for cp = 1:p
    Hs(:,cp) = ffunc(X(s,:)', [zeros(n,cp-1), ones(n,1), zeros(n,p-cp)]);
  end;
  for ck = 1:n
    for cp = 1:p
      for cn = 1:n
        Xs(ck,(cp-1)*n+cn) = Cbrev(ck,cn)*Hs(cn,cp);
      end;
    end;
  end;
  X2 = X2 + Xs'*Xs;
end;

D = [(eye(n*p) - A'*inv(A*A')*A)*X2; A];
bd = [zeros(n*p,1); b];
min_svd = min(abs(svd(D)));
if min_svd < 1e-10
  error('The constraints are probably splittable');
end;
wa = inv(D'*D)*D'*bd;

w = reshape(wa,n,p);
w = w / norm(w(:,end));

%=======================================================%
function Ab = fullrank(A)

Ab = [];
for cnt = 1:size(A,1)
  Ab = [Ab; A(cnt,:)];
  if det(Ab*Ab') < 1e-10
    Ab = Ab(1:end-1,:);
  end;
end;
