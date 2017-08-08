% Conjugate Project Gradient with reorthogonalisation:

function [lambda,alpha]=CG3(F1,Gi2,d,e)

% F1=B1*K1^-1*B1'+B2*KS*B2';
% Gi2=B2*R2;
% d=B2*KS*f2+B1*K1^-1*f1 ;
% e=-R2'*f2;

P=eye(size(Gi2,1))-Gi2*(Gi2'*Gi2)^-1*Gi2';
x=(Gi2*((Gi2'*Gi2)^-1))*e;
k=1;
r=d-F1*x;
w=P*r;
p=w;
n=(length(e)+length(d));
summe=0;
eps=1;

while   eps>= 10e-15 && k < n
    
    g=dot(w,w)/dot(p,F1*p);
    
    x=x+g*p;
    
    r=r-g*F1*p;
    
    w=P*r;
    
    summe=summe+((w'*F1*p)/(p'*F1*p))*p;
    
    p=w-summe;
    
    eps=norm(g*p);
    k=k+1;
end
alpha=(Gi2*((Gi2'*Gi2)^-1))'*-r;
lambda=x;
end
