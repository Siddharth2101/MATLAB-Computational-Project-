%{
basically, this script finds roots (where a function crosses zero) using 2 diff methods:

1. bisection: for equations like f(x)=0, it takes an interval [a,b] where the curve crosses the x-axis (meaning the ends have opposite signs). it chops the interval in half, checks which half the root is in, and repeats. keeps zooming in until the interval is smaller than the tol. slow but guaranteed to work.
2. fixed-point iteration: rewrites f(x)=0 as x=g(x). you take a starting guess, plug it into g(x) to get a new number, and keep plugging the new number back into g(x) over and over. some g(x) equations zero in on the root super fast, some are slow, and some just blow up and diverge entirely.
%}

clc
clear
close all

f1 = @(x) x-2.^(-x);
a1 = 0; 
b1 = 1; 
t1 = -8; 
to1 = 10^t1; 
m1 = 30; 
k1 = 1;
while k1 <= m1
    x1 = (a1+b1)/2;
    if f1(x1) == 0 || (b1-a1)/2 < to1
        break
    end
    if sign(f1(a1))*sign(f1(x1)) == -1
        b1 = x1;
    else 
        a1 = x1;
    end
    k1 = k1+1;
end
if k1 < m1+1
    fprintf('Tol: 10^%i, Approx: %.8f, Iters: %i\n', t1, x1, k1);
else
    fprintf('Max iters %i hit\n', m1);
end

f2 = @(x) x-2.^(-x);
a2 = 0; 
b2 = 1; 
t2 = -12; 
to2 = 10^t2; 
m2 = 30; 
k2 = 1;
while k2 <= m2
    x2 = (a2+b2)/2;
    if f2(x2) == 0 || (b2-a2)/2 < to2
        break
    end
    if sign(f2(a2))*sign(f2(x2)) == -1
        b2 = x2;
    else 
        a2 = x2;
    end
    k2 = k2+1;
end
if k2 < m2+1
    fprintf('Tol: 10^%i, Approx: %.12f, Iters: %i\n', t2, x2, k2);
else
    fprintf('Max iters %i hit\n', m2);
end

xp = 0:0.01:2.5;
figure
plot(xp, xp,'b','LineWidth',2)
hold on
plot(xp,2*sin(xp),'k','LineWidth',2)
axis([0 2.5 -0.5 2.5])
legend({'x','2sin(x)'},'Location','NorthWest')

f3 = @(x) x-2*sin(x);
a3 = 1; 
b3 = 3; 
t3 = -8; 
to3 = 10^t3; 
m3 = 100; 
xo = 999; 
k3 = 1;
while k3 <= m3
    xn = (a3+b3)/2;
    if f3(xn) == 0 || abs((xn-xo)/xn) < to3
        break
    end
    if sign(f3(a3))*sign(f3(xn)) == -1
        b3 = xn;
    else
        a3 = xn;
    end
    k3 = k3+1;
    xo = xn;
end
if k3 < m3+1
    fprintf('Tol: 10^%i, Approx: %.8f, Iters: %i\n', t3, xn, k3);
else
    fprintf('Max iters %i hit\n', m3);
end

x0 = 1; 
t4 = -10; 
m4 = 150;

g1 = @(x) (20*x+21/x^2)/21;
fpi(g1, x0, m4, t4);

g2 = @(x) x-(x^3-21)/(3*x^2);
fpi(g2, x0, m4, t4);

g3 = @(x) x-(x^4-21*x)/(x^2-21);
fpi(g3, x0, m4, t4);

g4 = @(x) (21/x)^(1/2);
fpi(g4, x0, m4, t4);

fprintf('From fastest to slowest: b, d, a (c does not converge). \n')

function fpi(g, x0, m, t)
    k = 1;
    fprintf('k: 0\tx0: %.10f\n', x0) 
    while  k <= m
        x = g(x0); 
        fprintf('k: %i \tx_%i: %.10f\t|x_%i-x_%i|: %.10f\n', k, k, x, k, k-1, abs(x-x0)) 
        if abs(x-x0) < 10^t
            break
        end
        k = k+1;
        x0 = x;
    end
    if k < m+1
        fprintf('Tol: 10^%i, Approx: %.' , t);
        fprintf(num2str(abs(t)+1));
        fprintf('f, Iters: %i\n\n', x, k);
    else
        fprintf('Max iters %i hit\n\n', m);
    end
end
