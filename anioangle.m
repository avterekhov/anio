function [w, ang] = anioangle(X, C)
% [w, ang] = anioangle(X, C)
%
% Determine the coefficients of a quadratic cost function such
% that the angle between the plain of the first two PCs of the
% experimental data X and the plain of the solutions of the
% optimization problem is minimal. The cost function has the form:
% J(x(1), ..., x(n)) = w(1,1) x(1)   + ... + w(n,1) x(n) +
%             + 1/2 * (w(1,2) x(1)^2 + ... + w(n,2) x(n)^2)
%
%
% The approximation method is described in:
%
% Terekhov, A. V., Pesin, Y. B., Niu, X., Latash, M. L., & Zatsiorsky,
% V. M. (2010). An analytical approach to the problem of inverse
% optimization with additive objective functions: an application to
% human prehension. Journal of mathematical biology, 61(3), 423-453.
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

[V, L] = eig(cov(X));
Cbrev = eye(size(C,2))-C'*inv(C*C')*C;
Cb = fullrank(Cbrev);

[k, J] = fminunc(@(k)(abs(subspace((Cb*diag(k/norm(k)))', ...
                                   V(:,1:size(Cb,1))))), ones(size(C,2),1), ...
                 optimset('Display', 'off', 'LargeScale', 'off', 'Algorithm', 'active-set', 'MaxFunEvals', 1e5));
k = k/norm(k);
d = -Cbrev*diag(k)*mean(X)';
w = [d, k];
ang = J;

%=======================================================%
function Ab = fullrank(A)

Ab = [];
for cnt = 1:size(A,1)
  Ab = [Ab; A(cnt,:)];
  if det(Ab*Ab') < 1e-10
    Ab = Ab(1:end-1,:);
  end;
end;
