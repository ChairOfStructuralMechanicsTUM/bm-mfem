
%%% divideModell
function  [substructure] = createSubstructure(femModel,dim,Boundary)

[elementArrayLeft, elementArrayRight]=divideElements(femModel.elementArray,dim,Boundary);

nodeArrayLeft  =  getNodesFromElements(elementArrayLeft);
nodeArrayRight =  getNodesFromElements(elementArrayRight);
substructure(1) = Substructure(nodeArrayLeft,elementArrayLeft);
substructure(2) = Substructure(nodeArrayRight,elementArrayRight);

% ELEMENTS

% divideElements function is placed in Core/Elements
%

% NODES

    function [ nodeArray] = getNodesFromElements( elementArray )
        tempNodeArray=[];
        for i=1:length(elementArray)
            
            Nodes=elementArray(i).getNodes;
            tempNodeArray=[tempNodeArray Nodes];
        end
        % delete duplicates and store unique nodes
        
        [~,sortedIds]=unique(tempNodeArray.getId); %[ids,index]
        nodeArray=tempNodeArray(sortedIds);
        
    end
end
