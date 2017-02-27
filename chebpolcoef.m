% Autor: Jeu Tonete
% Data: 28/12/2016
% 
% Funcao para o calculo dos coeficientes para interpolaçao usando polinomio de Chebyshev
% parametros: fname - funçao a ser interpolada, exempli de utilizaçao syms x; f(x) = sin(x);
%                 n - numero de coeficientes
%               rci - inicio da regiao de convergencia, padrao -1
%               rco - fim da regiao de convergencia, padrao 1
%                 m - grau do polinomio gerado
% Resposta:       c - coeficientes calculados
%               pol - polinomio resultante da interpolaçao
%


function [c,pol] = chebpolcoef(fname,n,rci,rco,m)
c = zeros(1,n);                         % pre alocaçao de memoria
syms w;                                 % definiçao para uso do simbolo
syms y;
syms x;

fnamey(y)=fname( (y*(rco-rci)+rco+rci)/2 );  % aplica o escalamento para regiao de convergencia escolhida
                                             % transformaçao x para y
for k = 0:n-1
    %Tk(w)=chebyshevT(k, w);            % armazena polinomio de chebychev de ordem k
    for j = 0:n
        thetaj=double((j+0.5)*pi/(n+1));
        %yj=cos(thetaj);                % ponto para a interpolaçao
        %c(k+1) = c(k+1) + fnamey(yj)*Tk(yj);
        c(k+1) = double(c(k+1)+ fnamey(cos(thetaj))*cos(k*thetaj));
                                        % calculo do coeficiente de chebyshev
    end
    c(k+1)=double(c(k+1)*2/(n+1));
end

Tm = chebyshevT(0:m-1, y);              % é calculado os polinomios de chebycheby
Tm = Tm.';                              % transpoe-se a matriz com os polinomios para a multiplicaçao com os coeficientes
pol(y) = c(1)*Tm(1)/2 + c(2:m)*Tm(2:m); % encontra-se o polinomio de interpolaçao resultante
pol(x) = pol((2*x-rco-rci)/(rco-rci));  % transformaçao inversa (y para x)
%pol(x) = vpa(pol(x));
pol(x) = expand(pol);                   % espande-se o polinomio para a forma mais simplificada
pol(x) = vpa(pol(x));                   % tira os numeros da forma racional n/n

t = rci:.01:rco;                        % plota graficos para o intervalo de convergencia
figure;
plot(t,fname(t),'b',t,pol(t),'r')       % plota a funçao original em azul, e a interpolaçao em vermelho
legend(char(fname),'chebyshev approximation')
figure;
plot(t,fname(t)-pol(t))                 % plota o erro entre a funçao original e a interpolada
legend('chebyshev approximation error')
% [c,p]=chebpolcoef(f,6,-1,1)
% t = rci:.01:rco;
% plot(t,fname(t),'b',t,pol(t),'r')
% plot(f(t)-p(t))

% >> f = @(x) 1./(1+25*x.^2);
% >> [c,x] = chebpolfit(f,11);
% >> t = -pi:.01:pi;
% >> plot(t,f(t),'b',t,chebpolval(c,t),'r',x,f(x),'ok')
% syms w
% Tn = chebyshevT([0, 1, 2, 3, 4], w)
% Tn = Tn.'
% pol(w) = c*Tn
% expand(pol(w))
% vpa(pol)
%  plot(t,f(t),'b',t,chebpolval(c,t),'r',t,pol(t),'-k')