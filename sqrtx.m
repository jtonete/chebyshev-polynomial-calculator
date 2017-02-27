syms x;

f(x)=sqrt(x);
[c,p]=chebpolcoef(f,30,0.5,1,5);

bits=15;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=1.454895*x-1.34491*x^2+1.106812*x^3-0.536499*x^4+0.1121216*x^5+0.2075806;

t = 0:.01:1;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 0.5:.01:1;
figure;
plot(t,abs(f(t)-p(t)),'r',t,abs(f(t)-prounded(t)),'b',t,abs(f(t)-faprox(t)),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')


