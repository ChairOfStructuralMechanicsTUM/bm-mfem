classdef IdTracker < handle
    %This class is used to track the maximum Ids of the nodes and of the
    %elements.
    
    properties (Access = private)
        nodeId
        elementId
    end
    
    methods
        %constructor
        function idTracker = IdTracker(nodeId, elementId)
            idTracker.nodeId = nodeId;
            idTracker.elementId = elementId;
        end
        %setter
        function setNodeId(idTracker, nodeId)
            idTracker.nodeId = nodeId;
        end
        function setElementId(idTracker, elementId)
            idTracker.elementId = elementId;
        end
        %getter
        function nodeId = getNodeId(idTracker)
            nodeId = idTracker.nodeId;
        end
        function elementId = getElementId(idTracker)
            elementId = idTracker.elementId;
        end
    end
    
end

