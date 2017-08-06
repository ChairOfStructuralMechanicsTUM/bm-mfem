
%%% divideModell
function  [substructure] = createSubstructure(femModel,dim,Boundary)

[elementArrayLeft, elementArrayRight,interfaceNodes]=divideElements(femModel.elementArray,dim,Boundary);

[nodeArrayLeft,interfaceNodeIds]  =  getNodesFromElements(elementArrayLeft,interfaceNodes);
substructure(1) = Substructure(nodeArrayLeft,elementArrayLeft);
substructure(1).interfaceNodes=interfaceNodeIds;

[nodeArrayRight,interfaceNodeIds] =  getNodesFromElements(elementArrayRight,interfaceNodes);
substructure(2) = Substructure(nodeArrayRight,elementArrayRight);
substructure(2).interfaceNodes=interfaceNodeIds;


% ELEMENTS

% divideElements function is placed in Core/Elements
%

% NODES

    function [ nodeArray,interfaceNodes] = getNodesFromElements( elementArray,interfaceNodes )
        tempNodeArray=[];
        
        for i=1:length(elementArray)           
            Nodes=elementArray(i).getNodes;
            tempNodeArray=[tempNodeArray Nodes];
        end
        % delete duplicates and store unique nodes
        
        [~,sortedIds]=unique(tempNodeArray.getId); %[ids,index]
        nodeArray=tempNodeArray(sortedIds);
        
  %      find Indices of interfaceNodes in nodeArray
                
        for i=1:length(interfaceNodes)
             interfaceNodes(i)=find(nodeArray.getId==interfaceNodes(i));
        end
        %Update location of InterfaceNodes:
        
        nodeArray=[nodeArray nodeArray(interfaceNodes)];
        nodeArray(interfaceNodes)=[];
        interfaceNodes=length(nodeArray)-length(interfaceNodes)+1:length(nodeArray);
    end
    
    %updateDofs(model.nodeArray,Substructure01,Substructure02,interfaceNodes)

        
    
end
