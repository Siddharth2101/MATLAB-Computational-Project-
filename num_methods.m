%{
file: num_methods_lab2.m

basically, here's how these two algorithms work:

1. robust quad roots (prob 1): finding roots for ax^2+bx+c but bulletproofed. first, it scales down the coeffs by dividing by their max val so we don't hit overflow/underflow with crazy big or small numbers. if a is 0, it just solves it as a linear eq. if the roots are real, it uses an alt version of the quadratic formula depending on the sign of 'b' to completely dodge catastrophic cancellation (which kills precision). returns a flag (0 for real, 1 for complex) plus the roots.

2. finite diff errors (prob 2): checking how step size (h) affects deriv approx for tan(x). calcs both forward difference and central difference, then compares them against the exact analytical deriv (sec^2(x)) to get the error. loops thru shrinking h values, finds where the error hits rock bottom before roundoff ruins it, and plots it all log-log to literally show the O(h) and O(h^2) slopes.
%}

t1 = rts(6, 5, -4);
if (t1(1) == 0)
    ty1 = 'real';
else
    ty1 = 'com';
end
fprintf("Solutions to 6x^2 + 5x - 4 = 0 are %s and are equal to %.5f and %.5f.\n \n", ty1, t1(2), t1(3));

t2 = rts(6*10^154, 5*10^154, -4*10^154);
if (t2(1) == 0)
    ty2 = 'real';
else
    ty2 = 'com';
end
fprintf("Solutions to (6*10^154)x^2 + (5*10^154)x - (4*10^154) = 0 are %s and are equal to %.5f and %.5f.\n \n", ty2, t2(2), t2(3));

t3 = rts(0,1,1);
if (t3(1) == 0)
    ty3 = 'real';
else
    ty3 = 'com';
end
fprintf("Solutions to x + 1 = 0 are %s and are equal to %.5f and %.5f.\n \n", ty3, t3(2), t3(3));

t4 = rts(1,-10^5,1);
if (t4(1) == 0)
    ty4 = 'real';
else
    ty4 = 'com';
end
fprintf("Solutions to x^2 - (10^5)x + 1 = 0 are %s and are equal to %.5f and %.5f.\n \n", ty4, t4(2), t4(3));

t5 = rts(1, -4, 3.999999);
if (t5(1) == 0)
    ty5 = 'real';
else
    ty5 = 'com';
end
fprintf("Solutions to x^2 - 4x + 3.999999 = 0 are %s and are equal to %.5f and %.5f.\n \n", ty5, t5(2), t5(3));

t6 = rts(10^(-155), 10^155, 10^155);
if (t6(1) == 0)
    ty6 = 'real';
else
    ty6 = 'com';
end
fprintf("Solutions to (10^-155)x^2 + (10^155)x + (10^155) = 0 are %s and are equal to %.5f and %.5f.\n \n", ty6, t6(2), t6(3));

t7 = rts(1, 1, 1);
if (t7(1) == 0)
    ty7 = 'real';
else
    ty7 = 'com';
end
fprintf("Solutions to x^2 + x + 1 = 0 are %s and have real and imaginary parts- %.5f and %.5f.\n \n", ty7, t7(2), t7(3));

x = 1;                      
k = 16;                     
plt(x,k);

function y = rts(a, b, c)
    s = max([abs(a), abs(b), abs(c)]);
    if s == 0
        error('All 0');
    end
    as = a / s;
    bs = b / s;
    cs = c / s;
    if as == 0
        rt = -cs / bs;
        y = [0, rt, rt];
        return;
    end
    dsq = bs^2 - 4*as*cs;
    if dsq < 0
        typ = 1;
        rp = -bs / (2 * as);
        ip = sqrt(-dsq) / (2 * as);
        y = [typ, rp, ip];
    else
        typ = 0;
        if bs > 0
            q = -bs-sqrt(dsq);
            x1 = 2*cs/q; 
            x2 = q/(2*as);
        else
            q = -bs+sqrt(dsq);
            x1 = q/(2*as); 
            x2 = 2*cs/q;
        end
        y = [typ, x1, x2];
    end
end

function y = fdt(x,h)
    y = (tan(x+h) - tan(x))/h;
end

function e = efd(x,h)
    e = abs(sec(x)^2 - fdt(x,h));
end

function y = cdt(x,h)
    y = (tan(x+h) - tan(x-h))/(2*h);
end

function e = ecd(x,h)
    e = abs(sec(x)^2 - cdt(x,h));
end

function [] = plt(x,k)
    h = 10.^-(0:k);            
    fd = zeros(size(h));
    cd = zeros(size(h));
    for i = 1:length(h)
        fd(i) = efd(x, h(i));
        cd(i) = ecd(x, h(i));
    end
    [~, fmi] = min(fd);
    [~, cmi] = min(cd);
    loglog(h, fd, 'r-o', 'LineWidth', 2, 'MarkerSize', 7);
    hold on
    loglog(h(2:7), 1e2*h(2:7).^1,'k --', 'LineWidth',2)
    loglog(h, cd, 'b--*', 'LineWidth', 2, 'MarkerSize', 7);
    loglog(h(3:6), 1e2*h(3:6).^2 - 10^(-5),'k -', 'LineWidth',2)
    loglog(h(fmi), fd(fmi), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y','LineWidth', 2);
    loglog(h(cmi), cd(cmi), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm','LineWidth', 2);
    xlabel('h');
    ylabel('Err');
    title('FD Err');
    legend('Fwd diff err', '$\mathcal{O}(h)$', 'Cent diff err', '$\mathcal{O}(h^2)$', 'Min FD err', 'Min CD err', 'Location', 'north', 'Interpreter', 'latex');
    set(gca,'Fontsize',10)
    xlim([h(end) h(1)])
    grid on
end
