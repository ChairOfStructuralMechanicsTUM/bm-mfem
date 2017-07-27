function [lambda, alpha]=CG2(F1,Gi2,d,e)
P=eye(4)-Gi2*(Gi2'*Gi2)^-1*Gi2';
x=(Gi2*((Gi2'*Gi2)^-1))*e;
r=d-F1*x;
w0=P*r;
p=w0;
%phi=0.5*x'*F1*x-d'*x;
b=0;
k=1;
n=(length(e)+length(d))+1; % usually convergence is reached in n steps
eps=1;
while  k < n  || eps > 10e-15 
    
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
    
      %  diff=abs(phi)-abs(0.5*x'*F1*x-d'*x)
      %  phi= 0.5*x'*F1*x-d'*x;
      % step(k)= g;
 
    eps=norm(g*p);    % => 
    k=k+1;
     

     % eps1=norm(d-[F1, -Gi2]*[x;alpha])  % => ||b-A*x||
end
 alpha=Gi2_inv'*-r;
lambda=x;


%loglog(step)
end