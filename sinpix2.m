syms x;
f(x)=sin(pi*x/2);

[c,p]=chebpolcoef(f,20,0,1,6);

bits=15;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=1.57035062*x+0.00508719*x^2-0.66666099*x^3+0.03610310*x^4  + 0.05512166*x^5;

t = -1:.01:1;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 0:.01:0.5;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')