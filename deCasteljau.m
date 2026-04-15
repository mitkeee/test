function y=deCasteljau(f,n,a,b,x)
% De Casteljauov algorithm, that evaluates the value of bernstein aprox.
% polynomial of degree n for a given function f on interval [a,b] 
% in a point x 

B=f(a+[0:n]/n*(b-a));
kon2=(x-a)/(b-a);
kon1=(1-kon2);
for r=1:n
    for i=0:n-r
        B(i+1)=kon1*B(i+1)+kon2*B(i+2);
    end
end
y=B(1);
end











