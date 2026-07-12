%{
basically, this script compares 3 iterative ways to solve Ax=b:

1. jacobi: splits matrix a into its diagonal and everything else. takes a guess for x, then uses that entire old guess to calc the new x all at once. it works but it's the slowest.
2. gs (gauss-seidel): splits matrix a into lower and upper triangles. unlike jacobi, as soon as it calculates a new piece of x, it uses it immediately for the next piece in the same loop. way faster.
3. sor (successive over-relaxation): takes gs and adds a cheat code (omega) to over-correct the guess, pushing it closer to the real answer even faster. the script tests a few omegas and then calcs the actual optimal one using the spectral radius.

the spectral radius is just the biggest absolute eigenvalue of the iteration matrix. the smaller it is, the faster the method converges to the answer.
%}

a = [ 10, -1, 2, 0; -1, 11, -1, 3; 2, -1, 10, -1; 0, 3, -1, 8];
b = [6; 25; -11; 15];
mit = 100;
tol = 10^-4;
x0 = [0; 0; 0; 0];

fprintf('\n -------Sol 1-------\n\n');

[ij, xj] = jac(a, b, mit, tol, x0);
fprintf('Jac:\nIters: %d\n', ij);
disp(xj);

[is, xs] = gsl(a, b, mit, tol, x0);
fprintf('\nGS:\nIters: %d\n', is);
disp(xs);

[mj, mg] = itr(a);
rj = max(abs(eig(mj)));
rg = max(abs(eig(mg)));
fprintf('\nRho Jac: %.6f\n', rj);
fprintf('Rho GS: %.6f\n', rg);

a2 = [ 4, -1, 0, 0, 0, 0; -1, 4, -1, 0, 0, 0; 0, -1, 4, 0, 0, 0; 0, 0, 0, 4, -1, 0; 0, 0, 0, -1, 4, -1; 0, 0, 0, 0, -1, 4];
b2 = [0; 5; 0; 6; -2; 6];
mi2 = 100;
to2 = 10^-8;
xo2 = zeros(6,1);

fprintf('\n -------Sol 2------\n\n');

[i2j, x2j] = jac(a2, b2, mi2, to2, xo2);
fprintf('Jac Iters: %d\n', i2j);

[i2s, x2s] = gsl(a2, b2, mi2, to2, xo2);
fprintf('GS Iters: %d\n', i2s);

omg = [1.01, 1.05, 1.1, 1.15];
xso = zeros(size(a2,1), length(omg));
its = zeros(1, length(omg));

for i = 1:length(omg)
    [it, x] = sor(a2, b2, omg(i), mi2, to2, xo2);
    xso(:,i) = x;
    its(i) = it;
end

fprintf('\nSOR Res:\n\n');
for i = 1:length(omg)
    fprintf('omg = %1.4f: Iters = %d\n', omg(i), its(i));
end

[m2j, ~] = itr(a2);
r2j = max(abs(eig(m2j)));
omo = 2/(1 + sqrt(1 - r2j^2));

[io, xo] = sor(a2, b2, omo, mi2, to2, xo2);
fprintf('\nOpt omg = %1.4f: Iters = %d\n', omo, io);
fprintf('Jac slowest, GS faster, SOR fastest.\n');

function [it, x] = jac(a, b, mit, tol, x0)
    n = diag(diag(a));
    p = n-a;
    it = 1;
    x = n\(p*x0+b);
    while it <= mit && norm(b - a*x, 2) / norm(b, 2) >= tol
        it = it+1;
        x0 = x;
        x = n\(p*x0+b);
    end
end

function [it, x] = gsl(a, b, mit, tol, x0)
    n = tril(a);
    p = n - a;
    it = 1;
    x = n\(p*x0 + b);
    while it <= mit && norm(b - a*x, 2) / norm(b, 2) >= tol
        it = it + 1;
        x0 = x;
        x = n\(p*x0 + b);
    end
end

function [it, x] = sor(a, b, om, mit, tol, x0)
    d = diag(diag(a));
    l = tril(a, -1);
    u = triu(a, 1);
    n = d + om*l;
    p = (1-om)*d - om*u;
    it = 1;
    x = n\(p*x0 + om*b);
    while it <= mit && norm(b - a*x, 2) / norm(b, 2) >= tol
        it = it + 1;
        x0 = x;
        x = n\(p*x0 + om*b);
    end
end

function [mj, mg] = itr(a)
    d = diag(diag(a));
    l = tril(a,-1);
    u = triu(a,1);
    mj = -d\(l+u);
    mg = -(d+l)\u;
end
