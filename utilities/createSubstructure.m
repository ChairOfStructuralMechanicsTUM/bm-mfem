
%%% divideModell
function  substructure = createSubstructure(femModel,dim,Boundary)

[elementArrayLeft, elementArrayRight]=divideElements(femModel.elementArray,dim,Boundary);

nodeArrayLeft  =  getNodesFromElements(elementArrayLeft);
nodeArrayRight =  getNodesFromElements(elementArrayRight);
substructure(1) = Substructure(nodeArrayLeft,elementArrayLeft);
substructure(2) = Substructure(nodeArrayRight,elementArrayRight);

% ELEMENTS
   
%         function [elementArrayLeft, elementArrayRight]= divideElements(elementArray,dim,Boundary)
%           elementArrayRight=[];
%           elementArrayLeft=[];
%           
%           for ii=1:length(elementArray)
%               
%               currentNodes = elementArray(ii).getNodes;
%               currentCoords = currentNodes.getCoords;           %[ x1 x2 y1 y2 z1 z2]
%               currentCoords = [currentCoords(2*dim-1),currentCoords(2*dim)];
%               
%               if currentCoords(2) >= currentCoords(1)    % get highest coordinat of the current element
%                   highCoord= currentCoords(2);
%                   lowCoord  = currentCoords(1);           % to compare element position                %
%                   lowCoordPosition=1;                   %lowCoordPosition is needed to overwrite the Node in Case of Boundary position
%               else
%                   highCoord = currentCoords(1);
%                   lowCoord  = currentCoords(2);
%                   lowCoordPosition=2;
%               end
%               
%               if highCoord <= Boundary  && lowCoord < Boundary         % if current element is left or on Boundary
%                   elementArrayLeft = [elementArrayLeft elementArray(ii)]; % add element to left  elementArray
%                   
%               elseif highCoord == Boundary && lowCoord == Boundary
%                                            
%                       elementArrayLeft = [elementArrayLeft copyElement(elementArray(ii))];
%                       elementArrayRight = [elementArrayRight copyElement(elementArray(ii))];
%                       % Set Cross-Section-Area
%                       setCrossSectionArea(elementArrayRight(length(elementArrayRight)),...
%                           0.5*getCrossSectionArea(elementArray(ii)));
%                       setCrossSectionArea(elementArrayLeft(length(elementArrayLeft)),...
%                           0.5*getCrossSectionArea(elementArray(ii)));
%                       
%                       
%                       %overwrite  both Nodes of Copied Element: old Node = current
%                       %Node ; new Node = copy of current Node
%                                          
%                       for j=1:2
%                           overwriteNode(elementArrayRight(length(elementArrayRight)),...
%                               currentNodes(j),copyElement(currentNodes(j)));
%                       end
%                       
%                   
%                   
%               elseif highCoord > Boundary && lowCoord ==Boundary
%                   elementArrayRight = [elementArrayRight elementArray(ii)];
%                     
%                 % overwrite Boundary node
%                      overwriteNode(elementArrayRight(length(elementArrayRight)),...
%                      currentNodes(lowCoordPosition),copyElement(currentNodes(lowCoordPosition)));
%                 
%               else
%                  elementArrayRight = [elementArrayRight elementArray(ii)]; 
%               end
%           end
%         end

% NODES

    function [ nodeArray] = getNodesFromElements( elementArray )
        tempNodeArray=[];
        for i=1:length(elementArray)
           
            Nodes=elementArray(i).getNodes;
            tempNodeArray=[tempNodeArray Nodes];
        end
        
       % tempNodeArray=[tempNodeArray InterfaceNodes];
        %[ids,index]
        [~,sortedIds]=unique(tempNodeArray.getId); %Deletion of duplicates
        nodeArray=tempNodeArray(sortedIds);     % only unique nodes are stored

        
        % Overwriting Id´s for Assembler ???
      
%          for ii=1:length(nodeArray)      
%               nodeArray(ii).id=ii;   
%          end
         
    end



end