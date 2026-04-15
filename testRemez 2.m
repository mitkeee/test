% test of the Remes algorithm
a=0;
b=2*pi;
tol=10^(-6);
fun=@(x) 36*exp(x/4).*sin(1.6*x);
maxk=11;

% approximating polynomial of degree 4:
% call the main Remez function
[P4,m,E0]=Remez(fun,a,b,4,tol,maxk);
disp('coefficients of approximating polynomial of degree 4:')
P4
%graph of the original function and the polinomial of the best uniform approximation:
xx=a:0.01:b;
plot(xx,fun(xx),'k--',xx,polyval(P4,xx),'g--')
legend('original function','polynomial of the best uniform approximation of degree 4')
xlabel('x')
ylabel('f(x)')
title('graph of the original function and the polinomial of the best uniform approximation of degree 4')
pause

% approximating polynomial of degree 6:
% call the main Remez function
[P6,m,E0]=Remez(fun,a,b,6,tol,maxk);
disp('coefficients of approximating polynomial of degree 6:')
P6
%graph of the original function and the polinomial of the best uniform approximation:
xx=a:0.01:b;
plot(xx,fun(xx),'k--',xx,polyval(P6,xx),'g--')
legend('original function','polynomial of the best uniform approximation of degree 6')
xlabel('x')
ylabel('f(x)')
title('graph of the original function and the polinomial of the best uniform approximation of degree 6')
pause