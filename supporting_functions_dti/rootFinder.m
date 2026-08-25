function [xi, maxiter, t0] = rootFinder(r, N)

[s] = size(r);
r = reshape(r, [prod(s), 1]);
% initial calculations ... 

lb = sqrt(2*N /corrterm(0, N) - 1);
maxiter = 500*ones(size(r));
epsi = 10e-12;


t0 = r - lb;
t1 = k(t0, N, r);

list = abs(t1-t0)>epsi;
while nnz(list)>0 
    t0(list) = t1(list);
    t1(list) = k(t0(list), N, r(list));
    maxiter(list) = maxiter(list)-1;
    
    t1(~isreal(t1))=0;
    if any(maxiter < 0)
        t0 = t1;
    end;
    
    list = abs(t1-t0)>epsi;
end; 

xi = corrterm(t1, N); xi(isnan(xi))=1;
maxiter =500 - maxiter;

xi = reshape(xi, s);
maxiter = reshape(maxiter, s);
end

function out = k(snr, N, r)
     snr(~isreal(snr))=0;
     snr(snr<0)=0;
    snr(isnan(snr))=100;
    term1 = snr;
    g_snr = g(snr, r, N);
    num = g_snr.*(g_snr - snr);
    try
        [K, dK] =kummer(-1/2,N,double(-(snr(:).^2)/2), 'fast');
    catch
        [K, dK] =kummer(-1/2,N,double(-(snr(:).^2)/2), 'slow');
    end
        %beta=sqrt(pi/2)*(doublefactorial(2*N-1)/((2^(N-1))*factorial(N-1)));
    beta = sqrt(2) * gamma(N+0.5)./ gamma(N);
    denum = snr.*(1+r.^2).*(1+(beta^2*K.*dK))-g_snr;
    out = term1 - num./denum;
end


function theta = g(snr, r, N)
    theta = sqrt(corrterm(snr, N).*(1+r.^2)-(2*N));
end


function xi = corrterm(snr, N)
    beta = sqrt(2) * gamma(N+0.5)./ gamma(N);
    %beta=sqrt(pi/2)*(doublefactorial(2*N-1)/((2^(N-1))*factorial(N-1)));
    snr(~isreal(snr))=0;
    snr(snr<0)=0;
    snr(isnan(snr))=100;
        K=kummer(-1/2,N,double(-snr(:).^2/2), 'slow');             
        xi = 2*N+snr.^2-beta^2*K.^2;
end

 