classdef DummyElement < Element
    %DUMMYELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        massMatrix
        stiffnessMatrix
        
    end
    
    methods
        
        function obj = DummyElement(id, nodeArray)
%             requiredPropertyNames = [];
            
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
        
        function check(obj)
           %CHECK override the element check function. nothing has to be
           %    checked for the dummy element.
        end
        
        function dofs = getDofList(obj)
            dofs([1 4]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 5]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3 6]) = obj.nodeArray.getDof('DISPLACEMENT_Z');
        end
    end
    
end

