n = 20;
h = 1;
x = 1;
deriv = cos(x);
hs = zeros(n,1);
errors = zeros(n,1);
for i = 1:n
    h = h/10;
    diffquo = (sin(x+h)-sin(x-h))/(2*h); % use central difference quotient
    error = abs(deriv - diffquo); % distance between actual derivative and the approximation
    hs(i) = h;
    errors(i) = error;
    fprintf('%g %g %g\n', h, diffquo, error)
end
loglog(hs, errors,'x')


function  = A()

end
