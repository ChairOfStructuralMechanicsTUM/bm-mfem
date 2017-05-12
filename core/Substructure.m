classdef Substructure < FemModel
    %Substructure The class of the individual substructure models
    %   This class is a substructure of the mesh in a FemModel. It can be 
    %   used to create substructes which are later used for FETI methods
    
    %??? FemModel or femModel in functions
    
    properties
    end
    
    methods
        %constructor
        function substructur = Substructure(nodeArray, elementArray, femModelParts)
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {nodeArray, elementArray};
            elseif nargin == 3
                super_args = {nodeArray, elementArray, femModelParts};
            end
            substructur@FemModel(super_args{:});
        end
        
        %function to divide FemModel into Substructures.
        %FemModel: structure to be divided
       
        %USE ELEMENTS INSTEAD OF NODES. NODES CAN BE GAINED FROM ELEMENTS
        %idsNodes: array of node ids along which structure should be
        %idsElements: array of element ids at the interface
        
        %divided. Need to be in order, smallest to biggest id.
        %substructure: array containing all the substructures
        %subs: index of current substructure
%         function substructure = divide(femModel, idsNodes, idsElements, substructure, subs)
%             
%             %-------------------2-------------------
%             %copy nodes along which Model is divided. 
%             cpsNodes = Substructure.copyNodes(idsNodes);
%             %copy elements which are at interface
%             cpsElements = Substructue.copyElements(idsElements);
%         
%             
%             %-------------------3-------------------
%             %find all nodes/elements that belong to one half or another in
%             %the substructure
%             %Substructure(1)
%             nodeArray = idsNodes;
%             elementArray = idsElements;
%             
%             %Substructure(2)
%             nodeArray = cpsNodes;
%             elementArray = cpsElements;
%             
%             
%             %-------------------4-------------------
%             %create actual substructures. call constructor for substruture
%             %to create it from now gained knowledge about nodes, elements
%             %and femParts included.
%             substructure.Substructure(nodeArray, elementArray, femModelParts); 
%            
%         end
        
       %copy nodes along which Model is divided. copies of nodes are
       %put in a new array called cps. this now contains all the
       %copies
       %ERROR: If there is an external force applied to the node, one needs
       %to also apply this load to the copied node (maybe only half of the
       %load?)
        function cpsNodes = copyNodes(idsNodes)
            cpsNodes = zeros(1,size(idsNodes,2));
            for node = 1:1:size(idsNodes,2)
                cpsNodes(node) = idsNodes(node).copyElement(idsNodes(node));
            end
        end
        
        %creates the elements at the interface from ids
        
        function cpsElements = copyElements(idsElements)
            cpsElements = zeros(size(idsElements,2));
            for element = 1:1:size(idsElements,2)
                cpsElements(node) = idsElements(element).copyElement(idsElements(node));
            end
        end
            
            
        
            
    end
        
        
end
    


