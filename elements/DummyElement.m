classdef DummyElement < Element
    %DUMMYELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        massMatrix
        dampingMatrix
        stiffnessMatrix
        dofArray
        restrictedDofs = Dof.empty
    end
    
    methods
        
        function obj = DummyElement(id, nodeArray)
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~ isa(nodeArray,'Node')
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, []};
            end
            
            % call the super class constructor
            obj@Element(super_args{:});
            obj.dofNames = [];
                        
        end
        
        function barycenter(obj)
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            stiffnessMatrix = obj.stiffnessMatrix;
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            massMatrix = obj.massMatrix;
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            dampingMatrix = obj.massMatrix;
        end
        
        function check(obj)
           %CHECK override the element check function. nothing has to be
           %    checked for the dummy element.
        end
        
        function getValuesVector(obj)
            
        end
        
        function update(obj)
            
        end
        
        function initialize(obj)
%             dofs = arrayfun(@(node) node.getDofArray, obj.nodeArray, 'UniformOutput', false);
%             obj.dofArray = [dofs{:}];
        end
        
        function dofs = getDofList(obj)
%             dofs = obj.dofArray(~obj.dofArray.isFixed);
            dofs = obj.dofArray;
%             dofs = [1:length(obj.dofArray)];
        end
        
        function setMatrices(obj, massMatrix, dampingMatrix, stiffnessMatrix)
            obj.massMatrix = massMatrix;
            obj.dampingMatrix = dampingMatrix;
            obj.stiffnessMatrix = stiffnessMatrix;
        end
        
        function setDofRestriction(obj, node, dofName)
        %SETDOFRESTRICTION if the dof is restricted to zero, it has to be
        %   removed from the dof list
%             obj.dofArray(obj.dofArray.getId==node.getDof(dofName).getId) = [];
%             node.fixDof(dofName);
        end
        
        function setDofOrder(obj, order)
%             dofs = Dof.empty;
            for ii=order
                obj.dofArray = [obj.dofArray obj.nodeArray(ii).getDofArray];
            end
%             obj.dofArray = dofs;
        end
    end
    
end

