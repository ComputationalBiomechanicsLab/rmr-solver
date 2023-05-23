function [A, jb, rows_rearranged, scaling_factors] = rref_custom(A, tol)
%RREF   Reduced row echelon form.
%   R = rref(A) produces the reduced row echelon form of A.
%
%   [R,jb] = rref(A) also returns a vector, jb, so that:
%       r = length(jb) is this algorithm's idea of the rank of A,
%       x(jb) are the bound variables in a linear system, Ax = b,
%       A(:,jb) is a basis for the range of A,
%       R(1:r,jb) is the r-by-r identity matrix.
%
%   [R,jb] = rref(A,tol) uses the given tolerance in the rank tests.
%
%   [R, jb, rows_rearranged] = rref_custom(A, tol)
%
%   Roundoff errors may cause this algorithm to compute a different
%   value for the rank than RANK, ORTH and NULL.
%
%   Class support for input A:
%      float: double, single
%
%   See also RANK, ORTH, NULL, QR, SVD.

%   Copyright 1984-2017 The MathWorks, Inc. 

[m,n] = size(A);

% Does it appear that elements of A are ratios of small integers?
[num, den] = rat(A);
rats = isequal(A, num./den);

% Compute the default tolerance if none was provided.
if (nargin < 2)
    tol = max(m,n)*eps(class(A))*norm(A,inf);
end

% Loop over the entire matrix.
i = 1;
j = 1;
jb = zeros(1,0);
rows_rearranged = 1:m;
scaling_factors = zeros(1,0);

while i <= m && j <= n
   % Find value and index of largest element in the remainder of column j.
   [p, k] = max(abs(A(i:m,j)));
   k = k+i-1;
   if p <= tol
      % The column is negligible, zero it out.
      A(i:m,j) = 0;
      j = j + 1;
   else
      % Remember column index
      jb = [jb j]; %#ok<AGROW>
      % Swap i-th and k-th rows.
      A([i k],j:n) = A([k i],j:n);
      rows_rearranged([i k]) = rows_rearranged([k i]);
      % Divide the pivot row by the pivot element.
      scaling_factors = [scaling_factors, A(i,j)];
      A(i,j:n) = A(i,j:n)./A(i,j);
      % Subtract multiples of the pivot row from all the other rows.
      for k = [1:i-1 i+1:m]
         A(k,j:n) = A(k,j:n) - A(k,j).*A(i,j:n);
      end
      i = i + 1;
      j = j + 1;
   end
end

% Return "rational" numbers if appropriate.
if rats
    [num, den] = rat(A);
    A = num./den;
end
