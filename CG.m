%Conjugate Project Gradient Method ( Feti_level_1)
function [lambda, alpha]=CG(F1,Gi2,d,e)

P=eye(size(Gi2,1))-Gi2*(Gi2'*Gi2)^-1*Gi2';
x=(Gi2*((Gi2'*Gi2)^-1))*e;
r=d-F1*x;
w0=P*r;
p=w0;
k=1;
n=(length(e)+length(d));
eps=1;
while  k < n  && eps > 10e-15
    
    g=dot(w0,w0)/dot(p,F1*p);
    
    x=x+g*p;
    
    r=r-g*F1*p;
    
    w=P*r;
    
    if norm(w)==0
        break
    end
    
    b=dot(w,w)/dot(w0,w0);
    
    w0=w;
    
    p=w+b*p;
    
    eps=norm(g*p) ;
    k=k+1;
end
alpha=(Gi2*((Gi2'*Gi2)^-1))'*-r;
lambda=x;
end