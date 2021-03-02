A = generateSPDmatrix(10);
b = rand(10,1);
xm = cgs(A,b)
[x, res, numItrs] = cg_iteration(A, b,  10^-17, length(b)*2) 

% CG on discrete Laplacian matrix
NValues = [16, 64, 256, 1024];
for N = NValues
    % G is an N+2 by N+2 full matrix with border entries set to zero
    G = numgrid('S',N+2);
    % A is the m by m five-point Laplacian sparse matrix, where m=N^2 
    A = delsq(G);
    m = N^2;
    b = ones(m,1);
    [x, res, numItrs] = cg_iteration(A, b, 10^-8, 1000);
    
    % graph it
    figure;
    semilogy(res);
    title(sprintf("N = %d", N));
end

% ||r_k|| does get reduced to 10^-8 in fewer than m = N^@ iterations

% functions
function [x, res, itrs] = cg_iteration(A, b, maxResidual, maxIterations)

    % set init vars
    itrs = 0;
    x = 0;
    r = b;
    p = r;
    res = [];
    
    while norm(r)>maxResidual && itrs < maxIterations % stopping conditions
        itrs = itrs + 1; % count number iterations

        Ap = A*p; % compute and reuse
        stepLength = (r'*r) / (p'*Ap);      % step length
        x = x + (stepLength * p);           % approx solution
        newR = r - (stepLength * Ap);       % residual
        beta = (newR'*newR) / (r'*r);       % improvement this step
        p = newR + (beta*p);                % search direction
        
        r = newR;
        res = [res norm(r)];
      
    end  

end

% taken from https://math.stackexchange.com/questions/357980/how-to-generate-random-symmetric-positive-definite-matrices-using-matlab
function A = generateSPDmatrix(n)
    % Generate a dense n x n symmetric, positive definite matrix
    A = rand(n,n); % generate a random n x n matrix
    % construct a symmetric matrix
    A = 0.5*(A+A'); 
    % since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
    % is symmetric positive definite, which can be ensured by adding nI
    A = A + n*eye(n);
end