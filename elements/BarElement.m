classdef (Abstract) BarElement < Element
    %BARELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        crossSectionArea
        length
    end
    
    methods
        % constructor
        function barElement = BarElement(id, material, crossSectionArea)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id; material};
            %added because needed for copy function
            elseif nargin == 2
                super_args = {id; material};
            end
            
            barElement@Element(super_args{:});
            
            %changed to 2 because of copy function
            if nargin > 2
                barElement.crossSectionArea = crossSectionArea;
            end
            
        end
        
        % getter functions
        function len = getLength(barElement)
            len = barElement.length;
        end
        
        function area = getCrossSectionArea(barElement)
            area = barElement.crossSectionArea;
        end
        
        % setter functions
        function setCrossSectionArea(barElement, area)
            barElement.crossSectionArea = area;
        end
        
        %%%Start NEW
        %function to half crosssection area at interface
        function halfCrossSectionArea(element)
            for ii = 1:length(element)
                element(ii).crossSectionArea = element(ii).crossSectionArea*0.5;
            end
        end
        
        %function preparing elements for copy and then sends them to copy
        function copiedElements = callToCopy(elementsToCopy, nodes01, nodes02, maxEleId, eleIntf)
            copiedElements = [];
            for ii = 1:length(elementsToCopy)
                nodePair = elementsToCopy(ii).getNodes;
                %new nodes need to be found for the new element, the nodes
                %are already copied before
                nodePair(1) = findCopy(nodePair(1), [nodes01 nodes02]); 
                nodePair(2) = findCopy(nodePair(2), [nodes01 nodes02]);

                
                if ismember(elementsToCopy(ii), eleIntf)    
                    copiedElements = [copiedElements copyElement(elementsToCopy(ii), maxEleId, nodePair)];
                     %increase element Id, so that no element has same id
                    maxEleId = maxEleId+1;
                else
                    id = elementsToCopy(ii).getId;
                    copiedElements = [copiedElements copyElement(elementsToCopy(ii), id, nodePair)];
                end
            end

        end
 
        %%%End NEW
        
        % member functions
        function update(barElement)
            barElement.length = computeLength(barElement.nodeArray(1).getCoords, ...
                    barElement.nodeArray(2).getCoords);
        end
        
        function c = barycenter(barElement)
            nodes = barElement.getNodes;
            c = (nodes(1).getCoords + nodes(2).getCoords) ./ 2;
        end
        
        function pl = draw(barElement)
            pl = line(barElement.nodeArray.getX, barElement.nodeArray.getY);
        end
        
        function pl = drawDeformed(barElement)
            pl = line(barElement.nodeArray.getX + barElement.nodeArray.getDofValue('DISPLACEMENT_X'), ...
                barElement.nodeArray.getY + barElement.nodeArray.getDofValue('DISPLACEMENT_Y'));
        end
        
        %function to sort elements in increasing order by their Id
        function elements = sortElements(elementsIn)
            elements = [];
            
            for ii = 1:length(elementsIn)
                mini = min(elementsIn.getId);
                ll = length(elementsIn);
                jj = 1;
                while jj <= ll
                    if elementsIn(jj).getId == mini
                        elements = [elements elementsIn(jj)];
                        elementsIn(jj) = [];
                        ll = ll-1;
                        jj = jj-1;
                    end
                    jj = jj+1;
                end
            end
        end
        
    end
    
    %%%Start NEW
    methods (Access = protected)
        %function that copies the elements
        function cp = copyElement(obj, maxId, cpNodes)
            %cp = copyElement@matlab.mixin.Copyable(obj);
            %create a basic BarElement
            cp = BarElement3d2n(maxId, obj.getMaterial);
            
            %instanciate node array like this, otherwise nodes get empty
            %degrees of freedeom. no BC is taken from original system
            cp.nodeArray = cpNodes;
            cp.crossSectionArea = obj.getCrossSectionArea;
            cp.length = computeLength(cpNodes(1).getCoords, cpNodes(2).getCoords);
        end
    end   
    %%%End NEW
end

