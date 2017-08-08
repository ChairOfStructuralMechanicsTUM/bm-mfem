
% Preconditioned Conjugate Project Gradtient
% (FETI_Level_1)
% Pi = ~ F1^-1  ( Preconditioner)
function [lambda,alpha]=PCG(Fi,Pi,Gi2,d,e)

P=eye(size(Gi2,1))-Gi2*(Gi2'*Gi2)^-1*Gi2';
x=(Gi2*((Gi2'*Gi2)^-1))*e;
r=d-Fi*x;
w=P*r;
%y=P*Pi*w;
p=w;

n=length(d)+length(e);
summe=0;
k=1;
eps=1;

while  eps >= 10e-15 && k < n
    
    g=(p'*w)/(p'*Fi*p);
        
    x=x+g*p;
    
    r=r-g*Fi*p;
     
    w=P*r;
    
    y=P*Pi*w;
    
    summe=summe+((y'*Fi*p)/(p'*Fi*p))*p;
    
    p=y-summe;
              
    k=k+1;
    eps=norm(p*g);
    
end

alpha=(Gi2*(Gi2'*Gi2)^-1)'*(Fi*x-d);
lambda=x;

end






