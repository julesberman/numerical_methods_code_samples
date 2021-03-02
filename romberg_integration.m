function [cur, k, maxm] = romberg(f, a, b)
    maxIter = 100;
    maxm = 0; % record the max value for m used
    % memoize functions for efficiency
    clearAllMemoizedCaches;
    memTkofm = memoize(@computeTkofm);
    memFforM = memoize(@evalFforM);
    
    startM = 4;
    prev = memTkofm(0,startM*2);
    cur = memTkofm(1,startM);
    k = 2;
    
    % main function
    while abs(cur-prev) > 10e-11 && k < maxIter
        prev = memTkofm(k-1,startM*2);
        cur = memTkofm(k,startM);
        k = k + 1;
    end
        
    % recursive functions    
    function Tkofm = computeTkofm(k,m)
        if k == 0
            % base case
            h = (b-a)/m;
            fm = memFforM(m);
            Tkofm = fm * h;
        else
            % richardson extrapolation
            fourtok = 4^k;
            Tkm1of2m = memTkofm(k-1,2*m);
            Tkm1ofm = memTkofm(k-1,m);
            Tkofm = (fourtok*Tkm1of2m - Tkm1ofm)/(fourtok - 1);
        end
    end

    function fm = evalFforM(m) % only works for powers of 2
        maxm = max(m, maxm);
        if m == 2
            % base case
            fm = ((f(a) + f(b)) / 2) + f((a+b)/2); % start with endpoints
        else
            h = (b-a)/m;
            fm = memFforM(m/2);
            xi = a + h;
            for ii = 1:2:m-1
                fm = fm + f(xi);
                xi = xi + (h*2);
            end
        end
    end    
end