syms x;

f(x)=log(x);
[c,p]=chebpolcoef(f,30,1,2,6);

bits=15;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=(0.9991150*(x-1)- 0.4899597*(x-1)^2+0.2856751*(x-1)^3-0.1330566*(x-1)^4+0.03137207*(x-1)^5);

t = 1:.01:2;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 1:.01:2;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')


