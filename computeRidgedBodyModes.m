
function [RidgedBodyModes,KS]= computeRidgedBodyModes(substructure)

StiffnessMatrix=SimpleAssembler(substructure).reducedStiffnessMatrix;
n=size(StiffnessMatrix,1);
L=zeros(n);
zeroPivots=[];     %Indices of  ZeroPivot Positions
nonZeroPivots=[];  % Indices of NonZeroPivot Positons
R=[];

% Cholesky factorization and Pivoting
    for j=1:n       
        sum1=0;       
        for k=1:j-1  
            sum1=sum1+(L(j,k))^2; 
        end
        diff=StiffnessMatrix(j,j)-sum1;
       
        if diff < 10^-5     % In case of Zero-Pivot
            zeroPivots(end+1)=j;
            % R=[R [StiffnessMatrix(1:j-1,j);1;zeros(n-j,1)]];
            
            
                        %  compute RBM  during factorization:
            %             y=L(1:j-1,1:j-1)\-r;
            %             rbm=L(1:j-1,1:j-1)'\y;
            %             rbm(j:n)=eye(n-j+1,1);
            %             RBM(:,length(ZZP))=rbm;

            % removing ZeroPivot variable          
            L(j,:)= 0;
            L(j,j)=1;
             continue
         
        else        % In case of Non Zero-Pivot
        L(j,j)=sqrt(StiffnessMatrix(j,j)-sum1);
        nonZeroPivots(end+1)=j;
        end
        
        for i =j+1:n
             sum2=0;
             for k=1:j-1
                 sum2=sum2+L(i,k)*(L(j,k));
             end
             L(i,j)=(StiffnessMatrix(i,j)-sum2)/L(j,j);
         end    
         
    end 
    
    % Reordering StiffnessMatrix
    Kpp=StiffnessMatrix(nonZeroPivots,nonZeroPivots);
         
         Kpr=StiffnessMatrix(nonZeroPivots,zeroPivots);
         Krr=StiffnessMatrix(zeroPivots,zeroPivots);
   
    
    
    % Compute Pseudo-Inverse
    KS=zeros(size(StiffnessMatrix));
    KS(1:size(Kpp,1),1:size(Kpp,2))= Kpp^-1;


      % Compute RBM   
    
    RidgedBodyModes=zeros(n,length(zeroPivots));
    
%     L=L(nonZeroPivots,nonZeroPivots);
%     n=size(L,1);  
%     for i = 1: length (zeroPivots)
%         R=StiffnessMatrix(nonZeroPivots,zeroPivots(i));
%         R(zeroPivots(i):n)=0;
%         y=L\-R;
%         rbm=L'\y;
%         RidgedBodyModes(1:zeroPivots(i),i)=rbm;
%         RidgedBodyModes(zeroPivots(i),i)=1;
%     end
    
    
    
    
    for i=1: length (zeroPivots)
        R=StiffnessMatrix(1:zeroPivots(i),zeroPivots(i));
        R(end)=1;
        if i ~=1
        R(zeroPivots(1:i-1))=0;    
        end
    y=L(1:zeroPivots(i),1:zeroPivots(i))\-R;%(1:zeroPivots(i),i);
    rbm=L(1:zeroPivots(i),1:zeroPivots(i))'\y;
    rbm(zeroPivots(i))=1;
    RidgedBodyModes(1:zeroPivots(i),i)= rbm;
    end
   
   
    % Kontrolle RBM 
   if norm(StiffnessMatrix*RidgedBodyModes)>10e-10 
       disp('StiffnessMatrix*RidgedBodyModes ~= 0 !!!')
   end
   if norm(Krr-Kpr'*Kpp^-1*Kpr) >101e-10
       disp('Krr-Kpr*Kpp^-1*Kpr ~= 0 !!! ')
   end
    
%   Kontrolle=K*KS*K - K   % if K = [Kpp Kpr;Kpr' Krr] 
end
