syms x;

f(x)=log(1+x);
[c,p]=chebpolcoef(f,30,0,1,6);

bits=15;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=(0.00001145 + 0.99916640*x- 0.48969909*x^2+0.28382318*x^3- 0.12995720*x^4+0.02980877*x^5);

t = 0:.01:1;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 0:.01:1;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')


