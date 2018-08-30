classdef FETI_DPSolver < handle
    
    
properties
end

methods
   
    function [lmd]=PCG(FIrr,FIrc,Kcc,Kccg,dr,fcg)
        %Initialisierinung
        lmd=zeros(length(dr),1);
        D=dr-FIrc*inv(Kccg)*fcg;
        A=(FIrr+FIrc*inv(Kccg)*FIrc.');
        n=length(lmd);
        
        di=FIrr*lmd;
        y=FIrc.'*lmd;
        x=inv(Kccg)*y;
        z=FIrc*x;
        di=di+z;

        r=D-di;
        d=r;
        x=lmd;
        while k<n || r>0.001 %Abbruchkriterium
            
          di=FIrr*x;
          y=FIrc.'*x;
          c=inv(Kccg)*y;
          z=FIrc*c;
          di=di+z;

          r=D-di;
            
          alpha= (r.'*r)/(d.'*A*d);
          x1=x+alpha*d;
          r1=r-alpha*A*d;
          beta=(r1.'*r1)/(r.'*r);
          d1=r1+beta*d;
          
          %Aufsetzten für nächsten Durchlauf
          x=x1;
          r=r1;
          d=d1;
          
          
          k=k+1;
  
        end
        lmd=x;
        
    end
    
    function [uc]=solveCornerDofs(Kccg,fcg,FIrc,lmd)
        uc=inv(Kccg)*(fcg+FIrc.'*lmd);
    end
    
    function[ur]=solveReminderDofs(Krr,gfr,Krc,Bc,uc,Br,lmd,hz,v)
        for i=1:hz
               for j=1:v
                   ur{j,i}=inv(Krr{j,i})*(gfr{j,i}-Krc{j,i}*Bc{j,i}*uc-Br{j,i}.'*lmd);
               end
        end
    end
    
 %alle Verschiebungen assemblen und zum plotten bereit machen:
    
end
    
    
    
    
    
    
    
end
    