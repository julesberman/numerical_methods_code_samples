%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION 5 ODE
% DE: θ''(t) = -sin(θ(t)) for 0 < t < T
% BC: θ(0)=alpha, θ(T)=beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Q5p1 %%%%%%%%%%%
T = 2*pi;
m = 32;
alpha = 0.7;
beta = 0.7;
% same solution as 2.4
u_guess = @(t)0.7*cos(t)+0.5*sin(t);
newtonIterate(T, m, u_guess, alpha, beta, 1);

% same solution as 2.5
u_guess = @(t)0.7 + sin(t/2);
newtonIterate(T, m, u_guess, alpha, beta, 1);

% % new solution
u_guess = @(t)2*cos(t/3)+sin(t).^2;
newtonIterate(T, m, u_guess, alpha, beta, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Q5p2 %%%%%%%%%%%
% same solution as 2.5 with larger T
m = 512;
T = 20;
u_guess = @(t)0.7 + 3*sin(t*1.4);
newtonIterate(T, m, u_guess, alpha, beta, 0);

% set larger Ts
T = 40;
u_guess = @(t)0.7 + 3*sin(t*1.4);
newtonIterate(T, m, u_guess, alpha, beta, 0); 

T = 120;
u_guess = @(t)0.7 + 3*sin(t*1.4);
newtonIterate(T, m, u_guess, alpha, beta, 0); 

T = 240;
u_guess = @(t)0.7 + 3*sin(t*1.4);
newtonIterate(T, m, u_guess, alpha, beta, 0); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function j = J(v, h) % compute Jacobian for some theta
    n = length(v);
    d = -2 + h^2*cos(v);
    sd = ones(n-1,1);
    j = diag(d)+diag(sd,1)+diag(sd,-1);
    j = j*1/h^2;
end

% compute G(theta)
function res = G(vk, h, alpha, beta)
    n = length(vk);
    res = zeros(n,1);
    hsq = h^2;

    for i = 1:n
        if i == 1
            vprev = alpha;
        else
            vprev = vk(i-1);
        end
        if i == n
            vnext = beta;
        else
            vnext = vk(i+1);
        end
        vi = vk(i);
        res(i) = ((vprev-(2*vi)+vnext)/hsq) + sin(vi);
    end
end

% take on step
function [vk1, Gt] = newtonStep(vk, h, alpha, beta)
    Gt = G(vk, h, alpha, beta);
    Jt = J(vk, h);
    update = Jt\(-Gt);
    vk1 = vk + update;
end

function [vk, normsGt, i] = newtonIterate(T, m, u, alpha, beta, make_plots)
    h = T/(m+1);
    % make t
    t = h:h:T-h;
    % transform t with guess of u function
    vk = u(t)';

    % set up
    i = 1;
    maxIter = 50;
    Gt = 1;
    vks = zeros(m,maxIter);
    normsGt = zeros(1,maxIter);
    vks(:,1) = vk; % record init
    % iterate with norm breaking condition
    while norm(Gt) > 10e-10 && i < maxIter    
        [vk, Gt] = newtonStep(vk, h, alpha, beta); 
        normsGt(i) = norm(Gt);
        vks(:,i+1) = vk;
        i = i + 1;
    end
    % truncate unused entries in arr
    normsGt = normsGt(1:i); 
    vks = vks(:,1:i); 
    
    % only plot last one
    if make_plots == 0
        
        figure;
        plot(t, vks(:,end), '.-', 'DisplayName', "Solution")
        ustr = func2str(u);
        ustr = ustr(5:end); %removes the '@(x) handle
        hold off; legend show, title(sprintf('Convergence of Newton Iterates, T=%0.3f, Init Guess: %s', T, ustr)), shg; 
        
        % plot norms of Gt
        figure; semilogy(normsGt), xlabel('iteration'), ylabel('norm(f(v)) (log scale)'), title('norm(f(v)) per iteration'), shg;
    end
    if make_plots == 1
        % do plots
        figure;
        for p = 1:i
           plot(t, vks(:,p), '.-', 'DisplayName', sprintf("itr %d", p-1))
           hold on
        end
        ustr = func2str(u);
        ustr = ustr(5:end); %removes the '@(x) handle
        hold off; legend show, title(sprintf('Convergence of Newton Iterates, T=%0.3f, Init Guess: %s', T, ustr)), shg; 
        
        % plot norms of Gt
        figure; semilogy(normsGt), xlabel('iteration'), ylabel('norm(f(v)) (log scale)'), title('norm(f(v)) per iteration'), shg;
    end
end