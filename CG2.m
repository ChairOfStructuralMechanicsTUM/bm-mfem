
% PCG nach FETI_Level_1

function [l,alpha]=CG2(Fi,Pi,Gi2,d,e)
P=eye(4)-Gi2*(Gi2'*Gi2)^-1*Gi2';
% k=0;
l=(Gi2*((Gi2'*Gi2)^-1))*e;
r0=d-Fi*l;
z0=Pi*r0;

% k=1;

b=0;
p=z0;  % p=s
p=P*p;
k=1;
eps=1;

while  k<10 && eps>10e-5
    
    gamma=(r0'*z0)/(p'*Fi*p);
    
    
    l=l+gamma*p;
    
    r=r0-gamma*Fi*p;
     
    z=Pi*r;
    
    b=dot(z,r)/dot(z0,r0);
       
    p=z+b*p;
    p=P*p;
    
        r0=r;
        z0=z;
        
    k=k+1;
    eps=norm(r);
    
end

alpha=(Gi2*(Gi2'*Gi2)^-1)'*(Fi*l-d);


end






