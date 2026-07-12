%{
basically, here's how these root-finding algorithms work:

1. secant method: it's a lot like newton's method but for when you don't actually have the derivative. it takes two initial guesses, draws a straight line (a secant) between them on the graph, and sees where it hits the x-axis to make the next guess. it repeats this till it zeros in on the root.
2. newton's method (scalar): uses the actual function and its derivative to slide down the tangent line straight to the root. usually, it's incredibly fast (quadratic convergence), but if the root is "flat" (multiplicity > 1, like in this problem), it slows down a lot and only converges linearly. the error plots literally show this drop in speed.
3. newton's method (systems): same exact vibe but for multiple variables at once. instead of dividing by a single derivative, it builds a jacobian matrix (a grid of all the partial derivatives) and solves a linear system at every single step to update all variables simultaneously.
%}

clc
clear all
close all

f1 = @(x) x.^2.*abs(sin(x))-4;
figure(1)
plot(0:0.01:4, f1(0:0.01:4),'b-','LineWidth',2)
xlabel('x'); ylabel('y')

p0 = 3.6; 
p1 = 3.7;
tle = -6; 
mit = 50; 
sct(p0, p1, f1, mit, tle)

p02 = 2.8; 
p12 = 2.9;
sct(p02, p12, f1, mit, tle) 

f2 = @(x) x.^2 - 2*x.*exp(-x) + exp(-2*x);
fpr = @(x) 2*(x - exp(-x)).*(1 + exp(-x));
x = 0:0.01:2;

figure(2)
plot(x, f2(x), 'b-', 'LineWidth', 2); hold on
plot(x, fpr(x), 'r-', 'LineWidth', 2)
legend('f(x)', 'fpr(x)')
xlabel('x'); ylabel('y')
grid on

tl2 = -10;          
mi2 = 50;
po2 = 1;
[pvs, n] = nts(f2, fpr, tl2, mi2, po2);
r = pvs(n);

etr = abs(pvs(1:n) - r);
kmx = length(etr) - 10;   
en = etr(1:kmx-1);
en1 = etr(2:kmx);

figure(3)
loglog(en, en1, 'o-','LineWidth',2); hold on
xlabel('en')
ylabel('en1')
grid on

i1 = floor(length(en)/2);
i2 = i1 + 5;
xrf = en(i1:i2);
c = en1(i1)/en(i1);
yrf = 1e+1*c * xrf;

figure(4)
loglog(xrf, yrf, 'k--','LineWidth',2)
legend('Dat','Slp=1','Location','best')

x0 = [0;0;0]; 
tl3 = 1e-6; 
mi3 = 50;
xs = nsy(@ffn, @jac, x0, tl3, mi3);
fprintf('\nFin:\n');
disp(xs)
ffn(xs)

function sct(p0, p1, f, mit, tle)
    tol = 10^tle;
    n = 2;   
    f0 = f(p0); 
    f1 = f(p1);
    fprintf('n: 1\tp%i: %.*f\t|f|: %.8f\t|err|: %.*f\n', n-1, abs(tle)+1, p1, abs(f1), abs(tle)+1, abs((p1-p0)/p0)) 
    while n <= mit
        p = p1-f1*(p1-p0)/(f1-f0);
        fprintf('n: %i\tp%i: %.*f\t|f|: %.8f\t|err|: %.*f\n', n, n, abs(tle)+1, p, abs(f(p)), abs(tle)+1, abs((p-p1)/p1)) 
        if abs((p-p1)/p1) < tol
             break
        end
        n = n+1;
        p0 = p1;
        f0 = f1;
        p1 = p;
        f1 = f(p);
    end
    if n < mit+1
        fprintf('Tol: 10^%i, App: %.10f, Iters: %i\n\n', tle, p, n);
    else
        fprintf('Max hit\n\n');
    end
end

function [pvs, n] = nts(f, fpr, tle, mit, p0)
    pvs = zeros(mit+1,1);
    pvs(1) = p0;
    fprintf('n: 0\tp0: %.*f\n', abs(tle)+1, p0)
    for n = 1:mit
        p = p0 - f(p0)/fpr(p0);
        pvs(n+1) = p;
        err = abs((p - p0)/p);
        fprintf('n: %d\tp%d: %.*f\t|err|: %.*f\n', n, n, abs(tle)+1, p, abs(tle)+1, err)
        if err < 10^tle
            break
        end
        p0 = p;
    end
    fprintf('\nTol: 10^%d, App: %.*f, Iters: %d\n\n', tle, abs(tle)+1, p, n);
end

function x = nsy(ffn, jac, x0, tol, mit)
    x = x0;
    err = 1.0;
    it = 0;
    fprintf('k\tx1\tx2\tx3\t||dx||\n');
    while err > tol && it < mit
        fx = ffn(x);
        jx = jac(x);
        dx = jx \ fx;
        x = x - dx;
        err = max(abs(dx));
        it = it + 1;
        fprintf('%d\t%.6f\t%.6f\t%.6f\t%.2e\n', it, x(1), x(2), x(3), err);
    end
end

function f = ffn(x)
    f = zeros(3,1);
    f(1) = 6*x(1) - 2*cos(x(2)*x(3)) - 1;
    f(2) = 9*x(2) + sqrt(x(1)^2 + sin(x(3)) + 1.06) + 0.9;
    f(3) = 60*x(3) + 3*exp(-x(1)*x(2)) + 10*pi - 3;
end

function j = jac(x)
    j = zeros(3,3);
    j(1,1) = 6;
    j(1,2) = 2*x(3)*sin(x(2)*x(3));
    j(1,3) = 2*x(2)*sin(x(2)*x(3));
    
    dnm = sqrt(x(1)^2 + sin(x(3)) + 1.06);
    j(2,1) = x(1)/dnm;
    j(2,2) = 9;
    j(2,3) = cos(x(3))/(2*dnm);
    
    j(3,1) = -3*x(2)*exp(-x(1)*x(2));
    j(3,2) = -3*x(1)*exp(-x(1)*x(2));
    j(3,3) = 60;
end
