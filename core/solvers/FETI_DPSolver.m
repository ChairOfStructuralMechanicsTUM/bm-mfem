classdef FETI_DPSolver < handle
    
    
properties
end

methods (Static)
   
    function [lmd]=PCG(FIrr,FIrc,Kccg,dr,fcg,lP)
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
        d=lP*r;
        x=lmd;
        k=1;
        while k<n | r>0.001 %Abbruchkriterium
            
          di=FIrr*x;
          y=FIrc.'*x;
          c=inv(Kccg)*y;
          z=FIrc*c;
          di=di+z;

          r=D-di;
            
          alpha= (r.'*lP*r)/(d.'*A*d);
          x1=x+alpha*d;
          r1=r-alpha*A*d;
          beta=(r1.'*lP*r1)/(r.'*lP*r);
          d1=lP*r1+beta*d;
          
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
    
    function[urem]=solveReminderDofs(Krr,gfr,Krc,Bc,uc,Br,lmd,v,hz)
        urem=cell(v,hz);
        for i=1:hz
               for j=1:v
                   urem{j,i}=inv(Krr{j,i})*(gfr{j,i}.'-Krc{j,i}*Bc{j,i}*uc-Br{j,i}.'*lmd);
               end
        end
    end
    
 %alle Verschiebungen assemblen und zum plotten bereit machen:
 function [ufinal]=getResultVector(femModel,uc,urem,ur,ur2,bcdof,v,hz)
     %Werte der dofs den richtigen globalen dof ids zuordnen
     urfinal=[];
           for i=1:hz
               for j=1:v
                   rDof=urem{j,i};
                   n=length(urfinal);
                   m=length(rDof);
                   urfinal(n+1:n+m)=rDof;  %alle dof verschiebungen der internen und boundry reminder dofs, inerface dofs doppelt enthalten
               end
           end
           
           %doppelte dofs an den interfaces eliminieren
           h=1;
           for l=1:length(ur2)
           helpVector=find(ur==ur2(l));
           if length(helpVector)==2
               dbnodes(h)=helpVector(1,2); %Einträge des globalen ur Nodes, die doppelt vorkommen
               h=h+1;
           end   
           end
           urfinal(dbnodes)=[];
           urfinal2=urfinal;

           %urfinal2=uniquetol(urfinal,0.0);  %toleranz ggf anpassen! Mögliche fehlerquelle!!!
           urIds=ur2.';
           %bcIds=unique(gbc,'stable') %Knotenids keine FGids!!!
           ubcIds=bcdof.getId;

           ufinal=zeros(length(femModel.getDofArray),1);
           ufinal(ubcIds)=uc;
           ufinal(urIds)=urfinal2;
           
           
        
           
 end
    
end
    
    
    
    
    
    
    
end
    