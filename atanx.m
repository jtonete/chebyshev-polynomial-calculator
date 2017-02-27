syms x;

f(x)=atan(x)/pi;
[c,p]=chebpolcoef(f,30,0,1,6);

bits=15;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=0.318253*x+0.003314*x^2-0.130908*x^3+0.068542*x^4-0.009159*x^5;

t = -1:.01:1;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = -1:.01:1;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')

t = 1:.1:10;
tinv = zeros(1,91);
for k = 1:91
    tinv(k)=1/t(k);
end
figure;
plot(t,f(t),'b',t,0.5-p(tinv),'r',t,0.5-prounded(tinv),'*y',t,0.5-faprox(tinv),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')
