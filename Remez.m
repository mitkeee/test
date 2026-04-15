function [P,m,E0]=Remez(f,a,b,n,tol,kmax)
%implementation of the Remese procedure (instructinos on e-classroom)
k=0;
endcondition=0;
E0=linspace(a,b,n+2);
while endcondition == 0 && (k<kmax)
    [P,m,endcondition,E1]=RemezStep (f,a,b,n,tol,E0)
    k=k+1;
    E0=E1,
end
end
%polynomials in matlab are given as vectors
function [P,m,endcondition,E1]=RemezStep(f,a,b,n,tol,E0)
%we use vandermond matrix
VV=vander(E0);
V=[(-1).^(0:n+1)',VV(:,2:end)];
sol=V\(f(E0))';
m=sol(1);
P=sol(2:end)';
%find extrema of the residual
xx=linspace(a,b,10^4+1);
[Mn,idx]=max(abs(f(xx)-polyval(P,xx)));
y=xx(idx);
if abs(Mn-abs(m))<tol
    endconndition=1;
    return
else
    endcondition=0;
    %change of points in set E0
    idx1=-1
    for j=1:n+1
        if (E0(j)-y)<=0&& (E0(j+1)-y>0), idx1=j;
        end
    end
    if y>=E0(end), idx1=n+2; end
    E1=E0;
    if idx1==-1 %y lies on [a,x_0]
        if sign(m)==sign(f(y)-polyval(P,y))
            E1(1)=y;
        else
            E1=[y,E0(1:(end-1))];
        end
    elseif idx1==n+2
     if sign((-1)^(n+1)*m)==sign(f(y)-polyval(P,y))
            E1(end)=y;
        else
            E1=[E0(2:(end),y)];
        end
    
 % is somewhere on the interval on the [x_{idx1},x_{idx1+1}]
elseif sign((-1)^(idx1+1)*m)==sign(f(y)-polyval(P,y))
        E1(idx1)=y;
    else
        E1=[idx1+1]==y;
    end
end
plot(xx,f(xx)-polyval(P,xx), 'k-')
hold on
plot(E0,0,'r*',y,0,'go',[y,y],[0,f(y)-polyval(P,y)],'g-')
plot(xx,zeros(1,length(xx)),'k--')
title('Graph of the residual')
pause
hold off
end