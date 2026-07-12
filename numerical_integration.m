%{
file: numerical_integration.m

basically, this script estimates the area under a curve (definite integral) using two methods:

1. composite trapezoidal rule: chops the area under the function into a bunch of trapezoids. it calculates the width of each chunk (h), adds up the endpoints once and all the interior points twice, and then multiplies the whole thing by h/2. it's simple but works decently well.
2. composite simpson's rule: instead of using straight lines, this fits little parabolas across sets of 3 points to build the area. it gives the endpoints a weight of 1, alternates weights of 4 and 2 for the interior points, and multiplies by h/3. it's usually way more accurate than the trapezoidal method for smooth curves.
%}

n = 6;
x = linspace(exp(1), exp(1)+2, n+1);
f = @(x) 1./(x.*log(x));

t = ctr(x, f(x));
s = csr(x, f(x));

iex = log(log(exp(1)+2));

fprintf('Trap Rule: %.6f\n', t)
fprintf('Simp Rule: %10.6f\n', s)
fprintf('Exct Val: %10.6f\n\n', iex)

fprintf('Abs Err (Trap): %.6e\n', abs(t - iex))
fprintf('Abs Err (Simp): %.6e\n', abs(s - iex))

function apx = ctr(x, f)
    h = (x(end)-x(1))/(length(x)-1); 
    apx = h/2*(f(1)+2*sum(f(2:end-1))+f(end)); 
end

function apx = csr(x, f)
    h = (x(end)-x(1))/(length(x)-1);
    apx = 1/3*h*(f(1)+4*sum(f(2:2:end-1))+2*sum(f(3:2:end-1))+f(end));
end
