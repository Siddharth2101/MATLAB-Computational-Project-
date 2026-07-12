%{
basically, this script tests two heavy-duty ways to estimate integrals (area under the curve) and compares their errors against the exact math answers:

1. clenshaw-curtis (cc): picks evaluation points by spacing them out using cosines so they bunch up near the edges (chebyshev nodes). it calculates the weights using trig formulas. it's super stable and great for most functions.
2. gauss-legendre (gl): the mathematically "perfect" way for polynomials. it uses the roots of legendre polynomials to pick points and finds the exact weights by calculating the eigenvalues of a specific tridiagonal matrix. it's designed to give zero error for polys up to degree 2n-1.

the script runs both methods on 4 diff types of functions, tracks how fast the error drops as you add more points (n), and plots it all on log scales so u can visually see which one converges better. 
%}

nmx = 50;
itr = [1/2; 2*(exp(-1) + sqrt(pi)*(erf(1)-1)); pi/2; 2/11];
lbl = {'|x|^3','exp(-x^{-2})','1/(1+x^2)','x^{10}'};

ecc = zeros(4, nmx);
ega = zeros(4, nmx);

for n = 1:nmx
    [xcc, wcc] = ccn(n);
    [xga, wga] = gau(n);
    
    fcc = tst(xcc);
    fga = tst(xga);
    
    for k = 1:4
        ecc(k,n) = abs(wcc * fcc(k,:)' - itr(k));
        ega(k,n) = abs(wga * fga(k,:)' - itr(k));
    end
end

plt(ecc, lbl, nmx, 'CC Errors');
plt(ega, lbl, nmx, 'GL Errors');

fprintf('\nDeg Exactness Test for GL:\n')
for n = 1:6
    [x,w] = gau(n);
    apx = w * (x.^10);
    err = abs(apx - 2/11);
    fprintf('N = %d, Err = %.2e\n', n, err);
end

function f = tst(x)
    f = zeros(4,length(x));
    f(1,:) = abs(x).^3;
    
    xtm = x;
    xtm(xtm==0) = eps;
    f(2,:) = exp(-xtm.^(-2));
    
    f(3,:) = 1./(1 + x.^2);
    f(4,:) = x.^10;
end

function plt(e, lbl, nmx, ttl)
    figure('Position',[100 100 900 900])
    sgtitle(ttl,'FontSize',18,'FontWeight','bold')
    for k = 1:4
        subplot(2,2,k)
        semilogy(1:nmx, e(k,:) + 1e-100,'o-','LineWidth',1.2)
        grid on
        xlim([1 nmx])
        ylim([1e-18 1e+1])
        xlabel('N','FontSize',14)
        ylabel('Err','FontSize',14)
        title(lbl{k},'FontSize',15,'FontWeight','bold')
    end
end

function [x,w] = ccn(n)
    tht = pi*(0:n)'/n;
    x = cos(tht);
    w = zeros(1,n+1);
    ii = 2:n;
    v = ones(n-1,1);
    if mod(n,2)==0
        w(1) = 1/(n^2 - 1);
        w(n+1) = w(1);
        for k = 1:(n/2 - 1)
            v = v - 2*cos(2*k*tht(ii)) / (4*k^2 - 1);
        end
        v = v - cos(n*tht(ii))/(n^2 - 1);
    else
        w(1) = 1/n^2;
        w(n+1) = w(1);
        for k = 1:((n-1)/2)
            v = v - 2*cos(2*k*tht(ii)) / (4*k^2 - 1);
        end
    end
    w(ii) = 2*v/n;
end

function [x,w] = gau(n)
    bta = .5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
    t = diag(bta,1) + diag(bta,-1);
    [v,d] = eig(t);
    x = diag(d);
    [x,ind] = sort(x);
    w = 2 * v(1,ind).^2;
end
