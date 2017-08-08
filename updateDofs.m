function updateDofs(nodeArray,Substructure01,Substructure02)
interfaceNodes=Substructure01.interfaceNodes;   
% fix  Dofs of copied Nodes on Boundary
% set CrossSectionArea on Boundary 

    for i =1:length(interfaceNodes)   
        
        currentNode=nodeArray(interfaceNodes(i));
        dofArray=currentNode.dofArray;
         
        index1=find(getId(Substructure01.nodeArray)  ==currentNode.id);
        index2=find(getId(Substructure02.nodeArray)  ==currentNode.id);

        for j=1:3
       % dofLoad(j)=dofArray(j).dofLoad
                
             if~isempty(dofArray(j).dofLoad)
           
            
            Substructure01.nodeArray(index1).dofArray(j).dofLoad=...
            0.5*dofArray(j).dofLoad  ;
            Substructure02.nodeArray(index2).dofArray(j).dofLoad=...
            dofArray(j).dofLoad;
            end
            if isFixed(dofArray(j))
            fix(Substructure01.nodeArray(index1).dofArray(j));
            fix(Substructure02.nodeArray(index2).dofArray(j));
            end
        end
    end
end
