%0.99998916 ? 0.99945630x+ 0.49556967x2-0.15375046x3 + 0.02553654x4
syms x;
f(x)=exp(-x);

[c,p]=chebpolcoef(f,25,0,1,5);

bits=16;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=vpa(0.99998916 - 0.99945630*x+ 0.49556967*x^2-0.15375046*x^3 + 0.02553654*x^4);

t = 0:.01:1;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 0:.01:1;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')