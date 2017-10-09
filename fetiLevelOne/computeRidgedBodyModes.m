
function [RidgedBodyModes,KS]= computeRidgedBodyModes(K)
n=size(K,1);
L=zeros(n);
zeroPivots=[];

% Cholesky factorization and Pivoting
for j=1:n
    sum1=0;
    for k=1:j-1
        sum1=sum1+(L(j,k))^2;
    end
    diff=K(j,j)-sum1;
    
    if diff < 10^-5     % In case of Zero-Pivot
        zeroPivots(end+1)=j;
        
        % removing ZeroPivot variable
        L(j,:)= 0;
        L(j,j)=1;
        continue
        
    else        % In case of Non Zero-Pivot
        L(j,j)=sqrt(K(j,j)-sum1);
        
    end
    
    for i =j+1:n
        sum2=0;
        for k=1:j-1
            sum2=sum2+L(i,k)*(L(j,k));
        end
        L(i,j)=(K(i,j)-sum2)/L(j,j);
    end
    
end

% Reordering StiffnessMatrix
%          Kpp=StiffnessMatrix(nonZeroPivots,nonZeroPivots);
%          Kpr=StiffnessMatrix(nonZeroPivots,zeroPivots);
%          Krr=StiffnessMatrix(zeroPivots,zeroPivots);


% Compute Pseudo-Inverse
KS=(L*L')^-1;
KS(zeroPivots,zeroPivots)= 0;

% Compute RBM
RidgedBodyModes=KS*(-K(:,zeroPivots));
for i=1:size(RidgedBodyModes,2)
    RidgedBodyModes(zeroPivots(i),i)=1;
end

% clean up "almost-zero" values in RBM:
for i=1:size(RidgedBodyModes,1)
    for j=1:size(RidgedBodyModes,2)
        if abs(RidgedBodyModes(i,j)) < 10e-15
            RidgedBodyModes(i,j)=0;
        end
    end
end
end
