m = 10;

% ground truth
A = diag(2*ones(1,m)) + diag(-1*ones(1,m-1),1) + diag(-1*ones(1,m-1),-1);
R = chol(A)

% should succeed and give same results as R above (in representation)
d = ones(1,m)*2;
s = ones(1,m-1)*-1;
[rd, rs] = tridia_cholesky(d,s)

% should throw error
s = ones(1,m-1)*-2;
[rd, rs] = tridia_cholesky(d,s)

%%%%%%%%%%%%%%%%%
% FLOP ANALYSIS %
%%%%%%%%%%%%%%%%%
% The algorithm uses about 5*(m-1) or O(m) floating point operations
% as there are 5 operations in the inter loop which runs m-1 times.
function [rd, rs] = tridia_cholesky(d,s)
    rd = zeros(size(d));
    rs = zeros(size(s));
    m = length(d);
    rd(1) = sqrt(d(1));
    for i = 1+1:m
        rs(i-1) = s(i-1)/rd(i-1);
        radicand = d(i) - rs(i-1)^2;
        if radicand < 0
            error("error: matrix is not postive definate")
        end
        rd(i) = sqrt(radicand);
    end
end

