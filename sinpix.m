syms x;
f(x)=sin(pi*x);

[c,p]=chebpolcoef(f,20,0,0.5,6);

bits=12;
prounded(x)=poly2sym(round(sym2poly(p*2^bits))/2^bits,x);
prounded(x)=vpa(prounded(x));

faprox(x)=3.140625*x+0.02026367*x^2-5.325196*x^3+0.5446778*x^4+1.800293*x^5;

t = -1:.01:1;
figure;
plot(t,f(t),'b',t,p(t),'r',t,prounded(t),'*y',t,faprox(t),'k')
legend('f(x)','chebyshev approx','rounded chebyshev','approx')

t = 0:.01:0.5;
figure;
plot(t,f(t)-p(t),'r',t,f(t)-prounded(t),'b',t,f(t)-faprox(t),'k')
legend('f(x)-chebyshev','f(x)-rounded chebyshev','f(x)-approx')

%points = [0:0.0625:0.5]
%plot(points,prounded(points)*2^12,'ob')