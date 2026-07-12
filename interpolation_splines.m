%{
file: interpolation_splines.m

basically, this script handles 3 ways to connect the dots (interpolation):

1. lagrange polys (prob 1): builds linear, quadratic, and cubic equations manually to fit a set of points for cos(log(x)). plots them all to show how higher degree polys get closer to the real curve, then spits out the absolute error at x=1.75.
2. divided differences (prob 2): uses newton's method to build a diff table and grab coeffs for an interpolating poly. tests it at specific points, plots it, and then adds a new data point to show how easily this method updates without starting over.
3. cubic splines (prob 3): fits smooth cubic curves piece-by-piece between data points (race times). builds a tridiag matrix to solve for the a,b,c,d coeffs. uses those to predict the time at the 0.75-mile mark and uses the derivative to find the speed at the start and finish lines.
%}

clc
close all
clear 

f = @(x) cos(log(x));
nod = [0.5 1 1.5 2];

p{1} = @(x) (x-nod(4))/(nod(3)-nod(4))*f(nod(3)) + (x-nod(3))/(nod(4)-nod(3))*f(nod(4));

p{2} = @(x) (x-nod(3)).*(x-nod(4))/((nod(2)-nod(3))*(nod(2)-nod(4)))*f(nod(2)) + (x-nod(2)).*(x-nod(4))/((nod(3)-nod(2))*(nod(3)-nod(4)))*f(nod(3)) + (x-nod(2)).*(x-nod(3))/((nod(4)-nod(2))*(nod(4)-nod(3)))*f(nod(4));

p{3} = @(x) (x-nod(2)).*(x-nod(3)).*(x-nod(4))/((nod(1)-nod(2))*(nod(1)-nod(3))*(nod(1)-nod(4)))*f(nod(1)) + (x-nod(1)).*(x-nod(3)).*(x-nod(4))/((nod(2)-nod(1))*(nod(2)-nod(3))*(nod(2)-nod(4)))*f(nod(2)) + (x-nod(1)).*(x-nod(2)).*(x-nod(4))/((nod(3)-nod(1))*(nod(3)-nod(2))*(nod(3)-nod(4)))*f(nod(3)) + (x-nod(1)).*(x-nod(2)).*(x-nod(3))/((nod(4)-nod(1))*(nod(4)-nod(2))*(nod(4)-nod(3)))*f(nod(4));

xx = 0.4:0.01:2.1;
apx = [1.75];

figure(1)
plot(xx, f(xx), 'k', 'LineWidth', 2)
hold on
for i = 1:3
    plot(xx, p{i}(xx), 'k--', 'LineWidth', 1)
    for j = 1:1
        fprintf('Appx at x=%.2f: %8.6f, Err: %.6f\n', apx(j), p{i}(apx(j)), abs(f(apx(j))-p{i}(apx(j))))
    end
end
plot(nod, f(nod), 'ro', 'MarkerSize', 8)
hold off
axis([0.4 2.1 0.5 1.1])

xd1 = [-0.1; 0.0; 0.2; 0.3];
fd1 = [17.3; 2.0; 5.19; 1.0];

cfs = div(xd1, fd1);
xv = [0.1 0.4];
ap1 = dvl(xd1, fd1, xv);
fprintf('Appx gen with %i pts\n', length(xd1))
for i = 1:length(xv)
    fprintf('Appx func val at %f is %.5f\n', xv(i), ap1(i))
end
fprintf('\n')

xv2 = -0.2:0.01:0.5;
fv1 = dvl(xd1, fd1, xv2);
figure(2)
plot(xd1,fd1,'or',xv2,fv1,'-k')

xd2 = [-0.1; 0.0; 0.2; 0.3; 0.05];
fd2 = [17.3; 2.0; 5.19; 1.0; 3.125];

ap2 = dvl(xd2, fd2, xv);
fprintf('Appx gen with %i pts\n', length(xd1))
for i = 1:length(xv)
    fprintf('Appx func val at %f is %.5f\n', xv(i), ap2(i))
end
fprintf('\n')

fv2 = dvl(xd2, fd2, xv2);
figure(3)
plot(xd1,fd1,'or',xv2,fv2)

x = [0 0.25 0.5 1 1.25];
a = [0 23.04 47.37 97.45 123.66];

[n, b, c, d] = csp(x, a);

x3 = 0.75;
i = 3;
tm = a(i)+b(i)*(x3-x(i))+c(i)*(x3-x(i)).^2+d(i)*(x3-x(i)).^3;
fprintf('Pred race tm: %i:%.2f\n\n', floor(tm/60), mod(tm, 60))

x4 = 0;
i = 1;
spd = 1/(b(i)+2*c(i)*(x4-x(i))+3*d(i)*(x4-x(i)).^2);
fprintf('Pred spd strt: %.2f mi/hr\n\n', spd*3600)

x5 = x(end);
i = n;
sp2 = 1/(b(i)+2*c(i)*(x5-x(i))+3*d(i)*(x5-x(i)).^2);
fprintf('Pred spd end: %.2f mi/hr\n\n', sp2*3600)

function [n, b, c, d] = csp(x, a)
    n = length(x)-1; 
    h = x(2:end)-x(1:end-1);
    am = zeros(n+1, n+1);
    am(2:end-1, 2:end-1) = diag(2*(h(1:end - 1) + h(2:end)))+diag(h(2:end-1), 1)+diag(h(2:end-1), -1);
    am(1,1) = 1;
    am(2,1) = h(1);
    am(end-1,end) = h(end);
    am(end,end) = 1;
    bm = zeros(n+1,1);
    bm(2:end-1) = 3./h(2:end).*(a(3:end)-a(2:end-1))-3./h(1:end-1).*(a(2:end-1)-a(1:end-2));
    c = (am\bm)';
    b = (a(2:end)-a(1:end-1))./h-h/3.*(2*c(1:end-1)+c(2:end));
    d = (c(2:end)-c(1:end-1))./(3*h);
end

function cfs = div(xd, fd)
    tb = zeros(length(xd), length(xd));
    tb(:, 1) = fd;
    for j = 2:length(xd)
        for i = j:length(xd)
            tb(i,j) = (tb(i,j-1)-tb(i-1,j-1))/(xd(i)-xd(i-j+1));
        end
    end
    cfs = diag(tb);
end

function apv = dvl(xd, fd, xv)
    tb = zeros(length(xd), length(xd));
    tb(:, 1) = fd;
    for j = 2:length(xd)
        for i = j:length(xd)
            tb(i,j) = (tb(i,j-1)-tb(i-1,j-1))/(xd(i)-xd(i-j+1));
        end
    end
    cfs = diag(tb);
    apv = zeros(1,length(xv));
    for i = 1:length(xv)
        x = xv(i);
        xdf = x-xd;
        prd = ones(1, length(cfs));
        for j = 2:length(cfs) 
            prd(j) = prod(xdf(1:j-1));
        end
        apv(1,i) = prd*cfs;
    end
end
