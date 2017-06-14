
function [RidgedBodyModes]= computeRidgedBodyModes(K)

K=K.reducedStiffnessMatrix;
n=size(K,1);
L=zeros(n);
zeroPivots=[];    %Indices of  ZeroPivot Positions
nonZeroPivots=[]; % Indices of NonZeroPivot Positons
R=[];

% Cholesky factorization and Pivoting
    for j=1:n       
        sum1=0;       
        for k=1:j-1  
            sum1=sum1+(L(j,k))^2; 
        end
        diff=K(j,j)-sum1;
       
        if diff < 10^-5     % In case of Zero-Pivot
            zeroPivots(end+1)=j;
            R=[R [K(1:j-1,j);1;zeros(n-j+1,1)]];
           
            
            %  compute RBM  during factorization:
%             y=L(1:j-1,1:j-1)\-r;

%             rbm=L(1:j-1,1:j-1)'\y;
%             rbm(j:n)=eye(n-j+1,1);
%             RBM(:,length(ZZP))=rbm;

            % removing ZeroPivot variable          
            L(j,:)=0;
            L(j,j)=1;
            continue
         
        else        % In case of Non Zero-Pivot
        L(j,j)=sqrt(K(j,j)-sum1);
        nonZeroPivots(end+1)=j;
        end
        
        for i =j+1:n
             sum2=0;
             for k=1:j-1
                 sum2=sum2+L(i,k)*(L(j,k));
             end
             L(i,j)=(K(i,j)-sum2)/L(j,j);
         end    
         
    end 
    
    % Reordering Stiffness 
    Kpp=K(nonZeroPivots,nonZeroPivots);
    Kpr=K(nonZeroPivots,zeroPivots);
    Krr=K(zeroPivots,zeroPivots);
  % Kontrolle=Krr-Kpr'*Kpp^-1*Kpr  
    
  %  Lreduced=L(NZZP,NZZP);  % L without ZeroPivot variable  
    
    % Compute Pseudo-Inverse
    KS=zeros(size(K));
    KS(1:size(Kpp,1),1:size(Kpp,2))= Kpp^-1;
    % Kontrolle=K*KS*K - K   % Wenn K in Form [Kpp Kpr;Kpr' Krr] vorliegt!
      
 
      % Compute RBM   
    RidgedBodyModes=zeros(n,length(zeroPivots));
    
    for i=1: length (zeroPivots)
    y=L(1:zeroPivots(i)-1,1:zeroPivots(i)-1)\-R(1:zeroPivots(i)-1,i);
    rbm=L(1:zeroPivots(i)-1,1:zeroPivots(i)-1)'\y;
    rbm(zeroPivots(i))=1;
    RidgedBodyModes(1:zeroPivots(i),i)= rbm;
    end
    
  % Kontrolle=K*RBM ;
    
end
