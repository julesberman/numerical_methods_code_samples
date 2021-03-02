clf;
n = 32;
A = diag(ones(n,1)*-1) + diag(ones(n-1,1),1) + diag(ones(n-2,1),2);

% from Def 2
t = 10;
grid = linspace(-t, t, 100);
L = length(grid);
Z2 = zeros(L);
X2 = zeros(L,1);
Y2 = zeros(L,1);

for j=1:L
    for k=1:L
        r = grid(j);
        im = grid(k);
        z = (r+im*1i);
        s = svd(A - eye(n)*z);
        X2(j) = r;
        Y2(k) = im;
        Z2(j,k) = s(end);
    end
end

contour(X2,Y2,Z2)
hold on
plotStyles = ["bx","rx","kx","b+","r+","k.","r.","b."];
% from Def 1
eLim = 8;
L = n;
X1 = zeros(L,1);
Y1 = zeros(L,1);
leg = {};
for j=1:eLim
    e = 10^-j;
    r = rand(n);
    dA = (r/norm(r))*e;
    A_pert = A + dA;
    eigs = eig(A_pert);
    
    for q = 1:length(eigs)
        eig_v = eigs(q);
        X1(q) = real(eig_v);
        Y1(q) = imag(eig_v);
     
    end
    plot(X1, Y1, plotStyles(j), 'LineStyle', 'none');

end
title('pseudospectrum of A')
xlabel('real part')
ylabel('img part')
axis equal
hold off
legend({'contours','10^-1','10^-2','10^-3','10^-4','10^-5','10^-6','10^-7','10^-8'})


% We see that the e-pseudospectrum for 10^-1 through 10^-8 is tightly
% compact around the origin. This means that the set of eigenvalues for
% matrices within a small pertubation of A grows largely relative to the
% pertubation. Thus small perturbations in A can result in very large 
% perturbations to the eigenvalues.