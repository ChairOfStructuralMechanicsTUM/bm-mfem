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
        %ids: array of node ids along which structure should be
        %divided. Need to be in order, smallest to biggest id.
        function substructure = divide(femModel, ids, substructure, subs)
            
            %if there are elements containing all the nodes, division is
            %possible otherwise an error occurs
            check = femModel.checkNodes(ids, femModel);
            if check == true
                fprintf('Choice of nodes for Substructe is fine.');
            else
                fprintf('FemModel cannot be substructured like this.');
                %substructuring is stopped
                return
            end
            
            %copy nodes along which Model is divided. 
            cps = femModel.copyNodes(ids)
        
            %find all nodes/elements that belong to one half or another in
            %the substructure
            
            
            %create actual substructures. call constructor for substruture
            %to create it from now gained knowledge about nodes, elements
            %and femParts included.
            substructure.Substructure(nodeArray, elementArray, femModelParts); 
           
        end 
        
        %function checks whether it is possible to divide geometry as
        %desired (all nodes in a row)
        function check = checkNodes(ids, femModel)
          
            for node = 1:1:size(ids,2)-1
                nodePair = [ids(node) ids(node+1)];
                %function elementArray in FemModel that gives out the whole
                %elementArray
                %??? ERROR
                check = ismember(nodePair, femModel.elementArray.barElement3d2n.nodeArray);
            end
             
            if check == 2*size(ids,2)
                check = true;
            else
                check = false;
            end           
        end
        
        
       %copy nodes along which Model is divided. dopies of nodes are
       %put in a new array called cps. this now contains all the
       %copies
       %ERROR: If there is an external force applied to the node, one needs
       %to also apply this load to the copied node (maybe only half of the
       %load?)
        function cps = copyNodes(ids)
            
            cps = zeros(1,size(ids,2));
            for node = 1:1:size(ids,2)
                cps(node) = ids(node).copyElement(ids(node));
            end
        end
        
            
    end
        
        
end
    


