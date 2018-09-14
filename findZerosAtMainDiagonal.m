
function findZerosAtMainDiagonal(M,K)

k=0;
for i = 1:length(M)
    if M(i,i) == 0
    fprintf('Main diagonal of M is zero at (%d,%d)',i)
    k=k+1;
    end
    if K(i,i) == 0
    fprintf('Main diagonal of M is zero at (%d,%d)',i)
    k=k+1;
    end    
end
if k == 0
    disp('No zeros at main diagonals')
end

