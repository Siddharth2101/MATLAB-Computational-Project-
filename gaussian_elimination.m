%{
file: gaussian_elimination.m

basically, this solves Ax=b using gaussian elimination, but uses pivoting so dividing by tiny numbers doesn't wreck the precision (cause float-pt math is sensitive). 

1. partial pivoting: merges a & b into one matrix. loops thru columns, finds the biggest absolute value in the current col, swaps that row to the top (the pivot), and zeros out everything underneath it. does this till it's an upper triangle, then just back-substitutes from bottom up to get the x vector.

2. scaled partial pivoting: same vibe, but first it finds the max absolute val in every single row to use as a 'scale'. when picking the best pivot, it divides the col values by their row's scale so it's looking for the largest *relative* number, not just the raw largest. swaps the rows AND the scales, eliminates, and back-subs.
%}

clc
a1 = [ 2  0  1  -1;
       6  3  2  -1;
       4  3 -2   3;
      -2 -6  2 -14];
 
b1 = [6; 15; 3; 12];
x1 = gpp(a1,b1);
fprintf('Sol vec\n');
fprintf('%8.4f\n', x1');

clc 
a2 = [     pi  -exp(1)   sqrt(2) -sqrt(3);
           pi   exp(1) -exp(1)^2      3/7;
      sqrt(5) -sqrt(6)         1 -sqrt(2);
         pi^3 exp(1)^2  -sqrt(7)     1/9];
b2 = [sqrt(11); 0; pi; sqrt(2)];
x2 = gsp(a2,b2);
fprintf('Sol vec\n');
fprintf('%8.4f\n', x2');

function x = gpp(a,b)
    a = [a b];
    n = size(a,1);
    for i = 1:n-1
        mid = find(abs(a(i:n, i)) == max(abs(a(i:n, i))));
        if mid(1) ~= 1 
            tmp = a(i, :);
            a(i, :) = a(mid(1)+i-1, :);
            a(mid(1)+i-1, :) = tmp;
            fprintf('Step %i: Swap %i & %i\n', i, i, mid(1)+i-1)
        else
            fprintf('Step %i: No swap\n', i)
        end
        for j = i+1:n
            a(j, i+1:end) = a(j, i+1:end)-a(j, i)/a(i, i)*a(i, i+1:end);
        end
        a(i+1:end, i) = zeros(n-i, 1);   
    end
    fprintf('\n')
    x = zeros(n, 1);
    x(n) = a(n, n+1)/a(n, n);
    for i = n-1:-1:1
        x(i) = (a(i, n+1)-dot(a(i, i+1:n), x(i+1:n)))/a(i, i);
    end
end

function x = gsp(a,b)
    scl = max(abs(a'))';
    a = [a b];
    n = size(a,1);
    for i = 1:n-1
        mid = find(abs(a(i:n, i)./scl(i:n)) == max(abs(a(i:n, i))./scl(i:n)));
        if mid(1) ~= 1
            tmr = a(i, :);
            a(i, :) = a(mid(1)+i-1, :);
            a(mid(1)+i-1, :) = tmr;
            tms = scl(i);
            scl(i) = scl(mid(1)+i-1);
            scl(mid(1)+i-1) = tms;
            fprintf('Step %i: Swap %i & %i\n', i, i, mid(1)+i-1)
        else
            fprintf('Step %i: No swap\n', i)
        end
        for j = i+1:n
            a(j, i+1:end) = a(j, i+1:end)-a(j, i)/a(i, i)*a(i, i+1:end);
        end
        a(i+1:end, i) = zeros(n-i, 1);   
    end
    fprintf('\n')
    x = zeros(n, 1);
    x(n) = a(n, n+1)/a(n, n);
    for i = n-1:-1:1
        x(i) = (a(i, n+1)-dot(a(i, i+1:n), x(i+1:n)))/a(i, i);
    end
end
