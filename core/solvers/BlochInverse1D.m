classdef BlochInverse1D < Solver
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        stiffnessMatrix
        massMatrix    
        
        leftDofs
        rightDofs
        leftNodes
        rightNodes     
    end
    
    methods
        function blochInverse1D=BlochInverse1D(femModel,leftNodes,rightNodes)
            if nargin > 0 
                blochInverse1D.femModel = femModel;
                blochInverse1D.assembler = SimpleAssembler(femModel);
                blochInverse1D.isInitialized = false;
            else
                error("Error (BlochInverse1D): no fem model defined!")
            end         
        end %blochInverse1D
        
        function solve(blochInverse1D, ~)
            if ~ blochInverse1D.isInitialized
                blochInverse1D.initialize();
            end
            
            
        end
       
        
        function initialize(blochInverse1D)
            if ~ blochInverse1D.femModel.isInitialized()
                blochInverse1D.femModel.initialize;
            end
            
            % assemble and reduce matrices
            blochInverse1D.massMatrix = blochInverse1D.assembler.assembleGlobalMassMatrix(blochInverse1D.femModel);
            blochInverse1D.stiffnessMatrix = blochInverse1D.assembler.assembleGlobalStiffnessMatrix(blochInverse1D.femModel);
            
            %%% suche nach weiteren Fehlern einfügen
            nodes = femModel.getAllNodes;
            
            a=2;

            while nodes(3,1) == nodes(3,a) %Compare y-Coordinates of first row
                a=a+1;
            end
            n=a-1; %length of the beam

            a=1;
            while nodes(2,1) == nodes(2,1+a*n) 
                %Compare x-Coordinates of first column

                a=a+1;
                if size(nodes,2) < 1+a*n
                    break
                end
            end
            m=a; %length (in y-direction) of the beam

            
       
            for i=1:m            
                x=1+n*(i-1);    %left node ID
                for j=1:n

                    if nodes(3,x)==nodes(3,x+j-1)
                    else
                        disp('y-Coordinates in the same row have to be equal')
                    end 
                end
            end

        end %initialize
        
        function [leftDofs,rightDofs]=getLeftRightDofs(NodeArray)
        
        end
        
      
        
      
        
    end %end methods
    
end