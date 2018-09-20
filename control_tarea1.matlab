%Resolver u(t) - Tj(t) - Tb(t) = 0
% O sea
% u(t) - JD2(theta)(t)  + bD(theta)(t) = 0
%CI = 0
 syms t
 rango = [0,3]
 y = dsolve('1 - D2y - Dy = 0','y(0)=0,Dy(0)=0')
 fplot(y,rango,'r:*')
 hold on
 %gráfica w(t)
 y2 = dsolve('1 - Dy - y = 0','y(0)=0')
 fplot(y2,rango,'b:+')
 %gráfica a(t)
 y3 = exp(-t)
 fplot(y3,rango,'g:s')
 hold off


