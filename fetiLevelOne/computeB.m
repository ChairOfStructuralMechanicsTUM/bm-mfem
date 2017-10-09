function [B]=computeB(substructure)

dofArray=substructure.getDofArray;
interfaceNodes=substructure.interfaceNodes;

%get only unfixed dof
fixedDofs = [];
for itDof = 1:length(dofArray)
    if dofArray(itDof).isFixed
        fixedDofs = [fixedDofs itDof];
    end
end
dofArray(fixedDofs) = [];

B=[];
for i =1: length(interfaceNodes)
    interfaceDofArray=substructure.nodeArray(interfaceNodes(i)).getDofArray;
    for j=1:3
        if ~isFixed(interfaceDofArray(j))
            B=[B;(dofArray==interfaceDofArray(j))];
        end
    end
end
end
