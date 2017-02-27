syms x;

f(x)=2*log10(x);
[c,p]=chebpolcoef(f,30,1,2,6);

bits=15;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
%pronuded(x)=vpa(prounded(x));

faprox(x)=(0.8678284*(x-1)-0.4255677*(x-1)^2+0.2481384*(x-1)^3-0.1155701*(x-1)^4+0.0272522*(x-1)^5);

t = 1:.01:2;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 1:.01:2;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')


