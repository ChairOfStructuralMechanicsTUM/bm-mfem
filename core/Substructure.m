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
            elseif nargin == 3
                super_args = {nodeArray; elementArray; femModelParts};
            end
            substructur@FemModel(super_args{:});
        end
        
        %function to divide FemModel into Substructures.
        %FemModel: structure to be divided
        %idsNodes: array of node ids along which structure should be
        %idsElements: array of element ids at the interface
        %divided. Need to be in order, smallest to biggest id.
        %substructure: array containing all the substructures
        %subs: index of current substructure
        function substructure = divide(femModel, idsNodes, idsElements, substructure, subs)
            
            %----------------1---------------------------------
            %not really necessary!
            
            %if there are elements containing all the nodes, division is
            %possible otherwise an error occurs
            check = Substructure.checkNodes(idsNodes, femModel);
            if check == true
                fprintf('Choice of nodes for Substructe is fine.');
            else
                fprintf('FemModel cannot be substructured like this.');
                %substructuring is stopped
                return
            end
            
            %-------------------2-------------------
            %copy nodes along which Model is divided. 
            cpsNodes = Substructure.copyNodes(idsNodes);
            %copy elements which are at interface
            cpsElements = Substructue.copyElements(idsElements);
        
            
            %-------------------3-------------------
            %find all nodes/elements that belong to one half or another in
            %the substructure
            %Substructure(1)
            nodeArray = idsNodes;
            elementArray = idsElements;
            
            %Substructure(2)
            nodeArray = cpsNodes;
            elementArray = cpsElements;
            
            
            %-------------------4-------------------
            %create actual substructures. call constructor for substruture
            %to create it from now gained knowledge about nodes, elements
            %and femParts included.
            substructure.Substructure(nodeArray, elementArray, femModelParts); 
           
        end 
        
        %function checks whether it is possible to divide geometry as
        %desired (all nodes in a row)
        function check = checkNodes(idsNodes, femModel)
          
            for node = 1:1:size(idsNodes,2)-1
                nodePair = [idsNodes(node) idsNodes(node+1)];
                %function elementArray in FemModel that gives out the whole
                %elementArray
                %??? ERROR
                check = ismember(nodePair, femModel.elementArray.barElement3d2n.nodeArray);
            end
             
            if check == 2*size(idsNodes,2)
                check = true;
            else
                check = false;
            end           
        end
        
        
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
    


