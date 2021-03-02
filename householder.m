%% 3) Householder
% If we were to compute (v_k)(v_k*) directly this would cost l^2 flops
% resulting in runtime of at least O(m^2*n^2)

% plot the performance
N = 1:2:100; % increments of n

% time arrays
Th = zeros(size(N));
Tq = zeros(size(N));

for i = 1 : length(N)
    % define A
    n = N(i);
    m = 10*n;
    A = rand(m,n);
    
    % time house
    tic;
    [W_h, R_h] = house(A);
    Th(i) = toc;
    
    % time qr
    tic;
    R_q = qr(A);
    Tq(i) = toc;
end

% plot house vs matlab-qr to see relationship
figure
hold on
xlabel('N')
ylabel('Time (s)')
title('QR Factorizations Times (house vs matlab-qr)')


plot(N,Th,'-o') % plot house
plot(N,Tq,'-x') % plot qr


legend('house','Matlab-QR')
hold off

% plot matlab-qr alone to more clearly see exponential runtime
figure
hold on
xlabel('N')
ylabel('Time (s)')
title('QR Factorizations Time (matlab-qr)')

plot(N,Tq,'-x') % plot qr

legend('Matlab-QR')
hold off

% house function
function [W,A] = house(A)
    dims = size(A);
    m = dims(1);
    n = dims(2);
    W = zeros(dims);
	for k = 1:n
        % set up x and e_1 of correct size
        x = A(k:end,k);
        e_1 = double(1:m-k+1 == 1)';
        
        % compute v_k
        v_k = sign(x(1))*norm(x)*e_1 + x;
        v_k = v_k / norm(v_k);
        
        % save v_k into W
        W(k:end,k) = v_k;
        
        % transform A to introduce zeros
        A(k:end,k:end) = A(k:end,k:end) - 2*v_k*(v_k'*A(k:end,k:end));
	end 
end

% formQ function
function [Q] = formQ(W)
    dims = size(W);
    m = dims(1);
    n = dims(2);
    Q = zeros(dims);
    for j = 1:n
        % create nth standard basis vector
        e_n = double(1:m == j)';
        
        % compute jth col of Q by multiplication with e_n
        for k = 1:n
            v_k = W(:,k);
            e_n(k:end) = e_n(k:end) - 2*v_k(k:end)*(v_k(k:end)'*e_n(k:end));      
        end
        
        % set jth col into Q
        Q(:,j) = e_n;
    end
end


