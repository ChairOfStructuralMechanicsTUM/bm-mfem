% Preconditioned CG (FETI_Level_1)
% Pi = ~ F1^-1  ( Preconditioner)

function [lambda,alpha]=PCG(Fi,Pi,Gi2,d,e)
% initializing
P=eye(size(Gi2,1))-Gi2*(Gi2'*Gi2)^-1*Gi2';
x=Gi2*(Gi2'*Gi2)^-1*e;
r=d-Fi*x;
w=P*r;
p=w;
n=length(d)+length(e);
summe=0;
k=1;
tol=1;

if w==0 % initial guess is already optimum
    alpha=(Gi2*(Gi2'*Gi2)^-1)'*(Fi*x-d);
    lambda=x;
else
    
    while   k < n && tol >= 10e-10
        
        g=(p'*w)/(p'*Fi*p);
        
        x=x+g*p;
        
        r=r-g*Fi*p;
        
        w=P*r;
        
        y=P*Pi*w;
        
        summe=summe+((y'*Fi*p)/(p'*Fi*p))*p;
        
        p=y-summe;
        
        k=k+1;
        
        tol=norm(p*g);
    end
    alpha=(Gi2*(Gi2'*Gi2)^-1)'*(Fi*x-d);
    lambda=x;
end
end






