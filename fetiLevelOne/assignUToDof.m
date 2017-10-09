function [U]=assignUToDof(substructure,u)
n=length(substructure.dofArray);
k=1;
U=zeros(length(u),2);
for i=1:n 
    if ~substructure.dofArray(i).fixed
        substructure.dofArray(i).value=u(k);
        
        U(k,:)=[substructure.dofArray(i).node.id, u(k)];
        k=k+1;
    end
end
