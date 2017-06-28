classdef EigensolverStrategy < Solver
    %EIGENSOLVERSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        stiffnessMatrix
        massMatrix
        
        eigenfrequencies
        modalMatrix
        spectralMatrix
    end
    
    methods
        
        function eigensolver = EigensolverStrategy(femModel)
           if nargin > 0
               eigensolver.femModel = femModel;
               eigensolver.assembler = SimpleAssembler(femModel);
               eigensolver.isInitialized = false;
           end
        end
        
        function solve(eigensolver)
            if ~ eigensolver.isInitialized
                eigensolver.initialize();
            end
            
            eigensolver.stiffnessMatrix = eigensolver.assembler.getStiffnessMatrix(eigensolver.assembler);
            eigensolver.massMatrix = eigensolver.assembler.assembleGlobalMassMatrix(eigensolver.femModel);
            
            %apply boundary conditions
            dofs = eigensolver.femModel.getDofArray();
            fixedDofs = [];
            for itDof = 1:length(dofs)
               if dofs(itDof).isFixed
                   fixedDofs = [fixedDofs itDof];
               end
            end
            eigensolver.stiffnessMatrix(fixedDofs,:) = [];
            eigensolver.stiffnessMatrix(:,fixedDofs) = [];
            eigensolver.massMatrix(fixedDofs,:) = [];
            eigensolver.massMatrix(:,fixedDofs) = [];
            
            % solve
            [eigensolver.modalMatrix, eigensolver.spectralMatrix] = ...
                eig(eigensolver.stiffnessMatrix, eigensolver.massMatrix);
            
            % get eigenfrequencies
            eigensolver.eigenfrequencies = diag(eigensolver.spectralMatrix);
            eigensolver.eigenfrequencies = sqrt(eigensolver.eigenfrequencies) ./ (2*pi);
        end
        
        function assignModeShapes(eigensolver)
            nModes = size(eigensolver.modalMatrix, 1);
            for itEv = 1:nModes
               eigensolver.assembler.assignResultsToDofs(eigensolver.femModel, eigensolver.modalMatrix.');
            end
        end
        
%         function modalSuperposition(eigensolver)
%             generalizedMassMatrix = eigensolver.modalMatrix.' ...
%                 * eigensolver.massMatrix ...
%                 * eigensolver.modalMatrix;
%             init_disp = [1 0];
%             
% %             generalizedStiffnessMatrix = eigensolver.modalMatrix.' ...
% %                 * eigensolver.stiffnessMatrix ...
% %                 * eigensolver.modalMatrix;
% 
%             xmu = eigensolver.modalMatrix * eigensolver.massMatrix * init_disp';
%             init = 1./diag(generalizedMassMatrix) .* xmu;
% 
% %             init = 1/diag(generalizedMassMatrix) * eigensolver.modalMatrix * eigensolver.massMatrix * init_disp;
%         end
        
        function initialize(eigensolver)
            eigensolver.femModel.initialize;
            eigensolver.isInitialized = true;
        end
        
        function ef = getEigenfrequencies(eigensolver)
            ef = eigensolver.eigenfrequencies;
        end
        
        function mm = getModalMatrix(eigensolver)
            mm = eigensolver.modalMatrix;
        end
        
    end
    
end

