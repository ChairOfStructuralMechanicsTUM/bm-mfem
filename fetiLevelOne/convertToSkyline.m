function [ sky, MAXA ] = convertToSkyline(A)

% % %Skyline:
 sky=[]; % store entries
 MAXA=[1 2]; % store indices
 n=length(A);
 for i=1:n
     firstNonZero=find(A(:,i),1);        % returns Indices of non-zero elements of 
     for j= i:-1:firstNonZero 
         sky(end+1) = A(i,j); 
     end
      MAXA(i+1)=MAXA(i)+(i-firstNonZero+1);
  end
  MAXA(end)=length(sky)+1;
end
con
