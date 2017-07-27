
%%% divideModell
function  [substructure,interfaceNodes] = createSubstructure(femModel,dim,Boundary)

[elementArrayLeft, elementArrayRight,interfaceNodes]=divideElements(femModel.elementArray,dim,Boundary);

[nodeArrayLeft, Id]  =  getNodesFromElements(elementArrayLeft,interfaceNodes);
substructure(1) = Substructure(nodeArrayLeft,elementArrayLeft);
substructure(1).interfaceNodes=Id;

[nodeArrayRight,Id] =  getNodesFromElements(elementArrayRight,interfaceNodes);
substructure(2) = Substructure(nodeArrayRight,elementArrayRight);
substructure(2).interfaceNodes=Id;


% ELEMENTS

% divideElements function is placed in Core/Elements
%

% NODES

    function [ nodeArray,Id] = getNodesFromElements( elementArray,interfaceNodes )
        tempNodeArray=[];
        for i=1:length(elementArray)
            
            Nodes=elementArray(i).getNodes;
            tempNodeArray=[tempNodeArray Nodes];
        end
        % delete duplicates and store unique nodes
        
        [~,sortedIds]=unique(tempNodeArray.getId); %[ids,index]
        nodeArray=tempNodeArray(sortedIds);
        % store Interface nodes at the end of nodeArray
        for i=1:length(interfaceNodes)
              Id(i)=find(nodeArray.getId==interfaceNodes(i));
% %             interfaceNode=nodeArray(Id);
% %             nodeArray(Id)=[];
% %             nodeArray=[nodeArray interfaceNode];
         end
         
        %Update location of InterfaceNodes:
        
%         b=length(nodeArray);
%         a=length(interfaceNodes);
%         a=b-a+1;
%         interfaceNodes=(a:b);
      
    end
    
    %updateDofs(model.nodeArray,Substructure01,Substructure02,interfaceNodes)

        
    
end
