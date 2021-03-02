%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3A
% Periodic Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interval [0,2pi]
L = 2*pi;

% create Ns to interate to grpah range
% use powers of 2
Ns = [];
for j = 1:7
    Ns = [Ns 2^j];
end
errs = [];
alpha = 1;
beta = 1;

for N = Ns
    % set up
    h = L/N;
    x = h*(0:N-1);
    
    % compute true f
    u_true = sin(cos(x));
    u_x_true = -sin(x).*cos(cos(x));
    u_xx_true = sin(x).^2.*(-sin(cos(x)))-cos(x).*cos(cos(x));
    f = u_xx_true + alpha*u_x_true - beta*u_true;

    % shift k coefficients for to match fft freqs
    k = fftshift((0:N-1)-(N/2));
    k((N/2)+1) = 0;
    
    % create factor in Fourier space and solve
    w = 1i*k;
    fact = (w.^2 + w.*alpha - 1*beta);
    f_hat = fft(f);
    u_hat = f_hat./fact;
    
    % inverse transform
    u_appx = real(ifft(u_hat));

    % record err
    errs = [errs norm(u_true-u_appx)];
    
end
semilogy(Ns, errs, '-x'); ylabel('err'); xlabel('N'); title("Periodic Boundary Conditions, err vs N pts");



%%%%%%%%%%%%%%%%%%%%%%
% PART 3C
% Sin Series Expansion
%%%%%%%%%%%%%%%%%%%%%%

% interval [0,pi]
L = pi;

% create Ns to interate to grpah range
% use powers of 2
Ns = [];
for j = 3:16
    Ns = [Ns 2^j];
end

errs = [];
alpha = 0;
beta = 1;

for N = Ns
    % set up
    h = L/N;
    x = h*(0:N-1);
    
    % compute true f
    u_true = sin(x).^2;
    u_xx_true = 2.*(cos(x).^2 - sin(x).^2);
    f = u_xx_true - beta*u_true;

    % build odd extension of f to force sine series
    % and in turn enforce boundry conditions 
    f_odd = f;
    f_odd(1) = 0;
    for i = 2:N
        f_odd(i) = -f(N-i+1);
    end
    f_odd = [f f_odd];
    
    % use N2 for odd extension
    N2 = N*2;
    
    % shift k coefficients for to match fft freqs
    k = fftshift((0:N2-1)-(N2/2));
    k((N2/2)+1) = 0;
    
    % create factor in Fourier space and solve
    w = 1i*k;
    fact = (w.^2 - 1*beta);
    f_hat = fft(f_odd);
    u_hat = f_hat./fact;
    
    % inverse transform
    u_appx = real(ifft(u_hat));
    u_appx = u_appx(1:N);
    
    % record err
    errs = [errs norm(u_true-u_appx)];
  
end
figure
loglog(Ns, errs, '-x'); ylabel('err'); xlabel('N'); title("Sine Series, err vs N pts");


