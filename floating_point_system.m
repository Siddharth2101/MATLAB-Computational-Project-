%{
file: floating_point_system.m

basically, here's how this models a float-pt sys:
1. gen elements: figures out every number a pc can rep with a given base, precision & exp range. pre-allocates mem for speed, builds a matrix for all mantissa combos, then loops thru exps, mantissas & leading digits. calcs the decimal, multiplies by base^exp, stores pos/neg versions, and sorts 'em.
2. ufl (underflow): smallest non-zero pos number before it just rounds to 0. it's just base^(lowest exp - 1).
3. ofl (overflow): absolute max val it can hold. max fraction * base^highest exp.
4. total & plot: uses basic p&c (2 signs * 1st digit options * mantissa options * total exps + 1 for zero) for the count. plot just drops the array on a 1D line so u can see the density.
%}

b = 2; 
t = 3; 
l = -1; 
u = 2;

els = gen(b,t,l,u);
disp(els);

n = num(b,t,l,u);
disp(n);

uf = ufl(b,l);
disp(uf);

of = ofl(b,t,u);
disp(of);

b2 = 10; 
t2 = 1; 
l2 = 1; 
u2 = 1;

el2 = gen(b2,t2,l2,u2);
disp(el2);

n2 = num(b2,t2,l2,u2);
disp(n2);

uf2 = ufl(b2,l2);
disp(uf2);

of2 = ofl(b2,t2,u2);
disp(of2);

plt(b,t,l,u);
plt(b2,t2,l2,u2);

function y = gen(b, t, l, u)
n = 2*(b-1)*(b^(t-1))*(u-l+1)+1;
y = zeros(n,1);
m = zeros(b^(t-1), t-1);
for k = 1:t-1
    for tm = 0:b^(t-k-1)-1
        for i = 0:b-1
            for j = 1:b^(k-1)
                m(tm*(b^k)+i*(b^(k-1))+j,k) = i;
            end
        end
    end
end
idx = 2;
for e = l:u
    for i = 1:b^(t-1)
        for j = 1:(b-1)
            s = 0;
            for k = 1:t-1
                s = s + ((1/b)^(k+1))*m(i,k);
            end
            s = s + (1/b)*j;
            s = s*(b^e);
            y(idx,1) = s;
            y(idx+1,1) = -1*s;
            idx = idx + 2;
        end
    end
end
y = sort(y);
end

function y = ufl(b, l)
y = b^(l-1);
end

function y = ofl(b, t, u)
y = (1 - b^(-t))*b^u;
end

function y = num(b, t, l, u)
y = 2*(b-1)*(b^(t-1))*(u-l+1)+1;
end

function [] = plt(b, t, l, u)
y = gen(b, t, l, u);
plot(y, 0, 'b.');
end
