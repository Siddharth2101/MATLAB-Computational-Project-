%{
file: matrix_conditioning.m

basically, this script tests how numerical rounding errors blow up when solving Ax=b for three types of matrices: hilbert, vandermonde, and tridiagonal. it sets a "true" x vector of all 1s, calculates b, and then forces matlab's backslash operator to solve for x again. then it compares the computed x to the real x to find the max error and grabs the matrix condition number (which measures how sensitive the matrix is to these errors). finally, it plots both error and condition number vs matrix size on loglog graphs. you'll see hilbert and vandermonde matrices get super ill-conditioned incredibly fast, while the tridiagonal one stays pretty stable. oh, and for vandermonde at n=8, the plot breaks because the error perfectly cancels out to exactly 0 due to float-pt quirks.
%}

phl(3,11);
pvn(4,17);
ptr(5,200);

function [e,k] = ehl(n)
    a = hilb(n);
    xt = ones(n,1);
    b = a*xt;
    xc = a\b;
    xe = abs(xc - xt);
    e = max(xe);
    k = cond(a);
end

function [] = phl(n1,n2)
    l = n2 - n1 + 1;
    nv = n1:n2;
    ev = zeros(1,l);
    kv = zeros(1,l);
    for i = 1:l
        [ev(i), kv(i)] = ehl(nv(i));
    end
    figure(1)
    subplot(1,2,1)
    loglog(nv, ev, 'r-o', 'LineWidth', 2, 'MarkerSize', 3);
    xlabel('dim'); ylabel('err');
    set(gca,'fontsize',14)
    subplot(1,2,2)
    loglog(nv, kv, 'b-o', 'LineWidth', 2, 'MarkerSize', 3);
    xlabel('dim'); ylabel('cond');
    set(gca,'fontsize',14)
    sgtitle('hilbert');
end

function [e,k] = evn(n)
    a = ones(n+1);
    for i = 0:n
        x = i/n;
        for j = 2:n+1
            a(i+1,j) = x^(j-1);
        end
    end
    k = cond(a);
    xt = ones(n+1,1);
    b = a*xt;
    xc = a\b;
    xe = abs(xc - xt);
    e = max(xe);
end

function [] = pvn(n1,n2)
    l = n2 - n1 + 1;
    nv = n1:n2;
    ev = zeros(1,l);
    kv = zeros(1,l);
    for i = 1:l
        [ev(i), kv(i)] = evn(nv(i));
    end
    figure(2)
    subplot(1,2,1)
    loglog(nv, ev, 'r-o', 'LineWidth', 2, 'MarkerSize', 3);
    xlabel('dim'); ylabel('err');
    set(gca,'fontsize',14)
    subplot(1,2,2)
    loglog(nv, kv, 'b-o', 'LineWidth', 2, 'MarkerSize', 3);
    xlabel('dim'); ylabel('cond');
    set(gca,'fontsize',14)
    sgtitle('vandermonde');
end

function [e,k] = etr(n)
    a = zeros(n);
    for i = 1:n
        a(i,i) = 2;
        if (i ~= 1)
            a(i,i-1) = -1;
        end
        if (i ~= n)
            a(i,i+1) = -1;
        end
    end
    k = cond(a);
    xt = ones(n,1);
    b = a*xt;
    xc = a\b;
    xe = abs(xc - xt);
    e = max(xe);
end

function [] = ptr(n1,n2)
    l = n2 - n1 + 1;
    nv = n1:n2;
    ev = zeros(1,l);
    kv = zeros(1,l);
    for i = 1:l
        [ev(i), kv(i)] = etr(nv(i));
    end
    figure(3)
    subplot(1,2,1)
    loglog(nv, ev, 'r-o', 'LineWidth', 2, 'MarkerSize', 3);
    xlabel('dim'); ylabel('err');
    set(gca,'fontsize',14)
    subplot(1,2,2)
    loglog(nv, kv, 'b-o', 'LineWidth', 2, 'MarkerSize', 3);
    xlabel('dim'); ylabel('cond');
    set(gca,'fontsize',14)
    sgtitle('tridiag');
end
