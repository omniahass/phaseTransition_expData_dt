function [mu, dmu] = chimean(nu, sig, N, selGr)
    flag = true;
    if nargout < 2
        flag = false;
        selGr = [false, false];
    end;
    selGr = selGr > 0;

    
    x = - nu.^2 ./ (2*sig.^2);
    beta = sqrt(2) * gamma(N+0.5)./ gamma(N);
    

    if ~flag
            [K] = kummer(-1/2, N, x, 'slow');
    else
            [K, dK] = kummer(-1/2, N, x, 'slow');
    end
        
        


    mu = beta .* sig .* K;

    
    if selGr(1)
        dxdnu = -nu./(sig.^2);
        dmudnu = beta .* sig.*dK .*dxdnu;
        
    else
        dmudnu = [];
    end
    if selGr(2)
        dxdsig = nu.^2 ./ sig.^3;
        dmudsig = beta*(K + sig.*dK.*dxdsig);
    else
        dmudsig = [];
    end;
    
    dmu = [dmudnu(:) dmudsig(:)];
end