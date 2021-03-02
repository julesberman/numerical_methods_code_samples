function [Amod,bmod]=gauss_el(A,b,display,pivot)
%
% On entry: A is a square matrix
%           b is a column vector of the same dimension
%           display is 1 if step-by-step display desired, 0 otherwise
% On exit: Amod and bmod are the result of applying Gauss elimination to
%          the system Ax=b. If pivot == 0, no row or column pivoting is done, so the
%          process is unstable and may break down even if A is nonsingular.
%          If pivot = 1, do partial pivoting (row interchanges)
%
n = length(b);
if size(A) ~= [n,n] | size(b) ~=[n,1],
   error('mismatched dimension')
end;
for k = 1:n-1,  % the kth step eliminates the kth column below the diagonal
              % element
    if pivot % partial pivoting (row interchanges)
    % find index q with largest entry in absolute value on or below the
    % diagonal in column k
        [maxval,q] = max(abs(A(k:n,k)));
        q = q + k-1; % want index between k and n
        A([k q],:) = A([q k],:); % swap rows k and q (ok if q == k)
        b([k q]) = b([q k]);
    end
    for i = k+1:n,  % We need to work on rows k+1 to n of the kth column
        mult = A(i,k)/A(k,k);  % multiplier which will give desired zero
        % A(i,:) = A(i,:) - mult*A(k,:); % subtract mult*(kth row) to ith row
        A(i,k+1:n) = A(i,k+1:n) - mult*A(k,k+1:n);  % don't need to compute the zeros
        A(i,k) = 0;  % so looks nice - or, could store mult here
        b(i) = b(i) - mult*b(k);       % apply same row operation to b
        if display
           A
           pause;
        end;
    end;
end;
Amod = A;
bmod = b;
