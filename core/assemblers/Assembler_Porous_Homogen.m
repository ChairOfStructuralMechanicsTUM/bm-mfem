classdef Assembler_Porous_Homogen < Assembler
    %Assembler_Porous_Homogen works as SIMPLEASSEMBLER, but is specified
    %for the application on mixed porous and homogen models
    
    properties
    end
    
    methods %test
        
        function assembling = Assembler_Porous_Homogen(femModel)
            if (nargin < 1)
                error('input model is missing');
            end
        end
    end
    
     methods (Static)
        
        function [K,Kred] = assembleGlobalStiffnessMatrix(femModel)
            mp1 = femModel.getModelPart('Homogen');
            mp2 = femModel.getModelPart('Porous');
            
            %create "missing" dofs for full matrix
            mp1.getNodes.addDof(["FRAME_DISPLACEMENT_X", "FRAME_DISPLACEMENT_Y", ...
            "FLUID_DISPLACEMENT_X", "FLUID_DISPLACEMENT_Y"]);
            mp2.getNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
            ndofs = length(femModel.getDofArray);
            
            node_array_mp1 = mp1.getElements.getNodeArray;
            node_array_mp2 = mp2.getElements.getNodeArray;
            node_array = femModel.getAllNodes;
            ids1 = mp1.getNodes().getId();
            ids2 = mp2.getNodes().getId();
            
            % matrix of homogenous part
            K_h = mp1.getElements.computeLocalStiffnessMatrix;
            
            K_h_full = zeros(ndofs);
            id_vector_homogen = zeros(size(K_h,1),1);
            for i_h=1:length(node_array_mp1)
                node = node_array_mp1(i_h);
                id_vector_homogen(2*i_h-1) = node.getDof('DISPLACEMENT_X').getId;
                id_vector_homogen(2*i_h) = node.getDof('DISPLACEMENT_Y').getId;
            end
            
            K_h_full(id_vector_homogen,id_vector_homogen)= K_h;
            
            % matrix of porous part
            K_p = mp2.getElements.computeLocalStiffnessMatrix;
            
            K_p_full = zeros(ndofs);
            id_vector_porous = zeros(size(K_p,1),1);
            for i_p=1:length(node_array_mp2)
                node = node_array_mp2(i_p);
                id_vector_porous(2*i_p-1) = node.getDof('FRAME_DISPLACEMENT_X').getId;
                id_vector_porous(2*i_p) = node.getDof('FRAME_DISPLACEMENT_Y').getId;
                id_vector_porous((size(K_p,1)/2)+2*i_p-1) = node.getDof('FLUID_DISPLACEMENT_X').getId;
                id_vector_porous((size(K_p,1)/2)+2*i_p) = node.getDof('FLUID_DISPLACEMENT_Y').getId;
            end
            
            K_p_full(id_vector_porous,id_vector_porous)= K_p;
    
            K_full=K_p_full+K_h_full;
            
            % Determining intersection (common nodes)
            coupling = intersect(ids1,ids2);

            % Applying coupling conditions
            delete = zeros(length(coupling),1);
            

            for iNode = 1:length(coupling)
                node_id = coupling(iNode);
                id_homogen_x = node_array(node_id).getDof('DISPLACEMENT_X').getId;
                id_homogen_y = node_array(node_id).getDof('DISPLACEMENT_Y').getId;
                id_frame_x = node_array(node_id).getDof('FRAME_DISPLACEMENT_X').getId;
                id_frame_y = node_array(node_id).getDof('FRAME_DISPLACEMENT_Y').getId;
                id_fluid_y = node_array(node_id).getDof('FLUID_DISPLACEMENT_Y').getId;

                K_full(id_homogen_x,:)=K_full(id_homogen_x,:)+K_full(id_frame_x,:);
                K_full(:,id_homogen_x)=K_full(:,id_homogen_x)+K_full(:,id_frame_x);
                K_full(id_homogen_y,:)=K_full(id_homogen_y,:)+K_full(id_frame_y,:);
                K_full(:,id_homogen_y)=K_full(:,id_homogen_y)+K_full(:,id_frame_y);
                K_full(id_homogen_y,:)=K_full(id_homogen_y,:)+K_full(id_fluid_y,:);
                K_full(:,id_homogen_y)=K_full(:,id_homogen_y)+K_full(:,id_fluid_y);

                delete([3*iNode-2,3*iNode-1,3*iNode]) = [id_frame_x,id_frame_y,id_fluid_y];
            end

            %reducing matrix
            for i_id1=1:length(ids1)
                if ~ismember(ids1(i_id1),coupling)
                    delete(end+1)=node_array_mp1(i_id1).getDof('FRAME_DISPLACEMENT_X').getId;
                    delete(end+1)=node_array_mp1(i_id1).getDof('FRAME_DISPLACEMENT_Y').getId;
                    delete(end+1)=node_array_mp1(i_id1).getDof('FLUID_DISPLACEMENT_X').getId;
                    delete(end+1)=node_array_mp1(i_id1).getDof('FLUID_DISPLACEMENT_Y').getId;
                end
            end

            for i_id2=1:length(ids2)
                if ~ismember(ids2(i_id2),coupling)
                    delete(end+1)=node_array_mp2(i_id2).getDof('DISPLACEMENT_X').getId;
                    delete(end+1)=node_array_mp2(i_id2).getDof('DISPLACEMENT_Y').getId;
                end
            end
            
            % apply fixed dofs
            [~, fixedDofs] = femModel.getDofConstraints;
                if ~ isempty(fixedDofs)
                    fixedDofIds = fixedDofs.getId();
                else
                    fixedDofIds = [];
                end

            K=K_full;
            K(delete,:) = [];
            K(:,delete) = [];
            
            Kred=K_full;
            Kred([delete' fixedDofIds],:) = [];
            Kred(:,[delete' fixedDofIds]) = [];  
        end
        
        function [M,Mred] = assembleGlobalMassMatrix(femModel)
            mp1 = femModel.getModelPart('Homogen');
            mp2 = femModel.getModelPart('Porous');
            node_array_mp1 = mp1.getElements.getNodeArray;
            node_array_mp2 = mp2.getElements.getNodeArray;
            node_array = femModel.getAllNodes;
            ids1 = mp1.getNodes().getId();
            ids2 = mp2.getNodes().getId();
            ndofs = length(femModel.getDofArray);
            
            % matrix of homogenous part
            M_h = mp1.getElements.computeLocalMassMatrix;
            
            M_h_full = zeros(ndofs);
            id_vector_homogen = zeros(size(M_h,1),1);
            for i_h=1:length(node_array_mp1)
                node = node_array_mp1(i_h);
                id_vector_homogen(2*i_h-1) = node.getDof('DISPLACEMENT_X').getId;
                id_vector_homogen(2*i_h) = node.getDof('DISPLACEMENT_Y').getId;
            end
            
            M_h_full(id_vector_homogen,id_vector_homogen)= M_h;
            
            % matrix of porous part
            M_p = mp2.getElements.computeLocalMassMatrix;
            
            M_p_full = zeros(ndofs);
            id_vector_porous = zeros(size(M_p,1),1);
            for i_p=1:length(node_array_mp2)
                node = node_array_mp2(i_p);
                id_vector_porous(2*i_p-1) = node.getDof('FRAME_DISPLACEMENT_X').getId;
                id_vector_porous(2*i_p) = node.getDof('FRAME_DISPLACEMENT_Y').getId;
                id_vector_porous((size(M_p,1)/2)+2*i_p-1) = node.getDof('FLUID_DISPLACEMENT_X').getId;
                id_vector_porous((size(M_p,1)/2)+2*i_p) = node.getDof('FLUID_DISPLACEMENT_Y').getId;
            end
            
            M_p_full(id_vector_porous,id_vector_porous)= M_p;
    
            M_full=M_p_full+M_h_full;
            
            % Determining intersection (common nodes)
            coupling = intersect(ids1,ids2);

            % Applying coupling conditions
            delete = zeros(length(coupling),1);
            

            for iNode = 1:length(coupling)
                node_id = coupling(iNode);
                id_homogen_x = node_array(node_id).getDof('DISPLACEMENT_X').getId;
                id_homogen_y = node_array(node_id).getDof('DISPLACEMENT_Y').getId;
                id_frame_x = node_array(node_id).getDof('FRAME_DISPLACEMENT_X').getId;
                id_frame_y = node_array(node_id).getDof('FRAME_DISPLACEMENT_Y').getId;
                id_fluid_y = node_array(node_id).getDof('FLUID_DISPLACEMENT_Y').getId;

                M_full(id_homogen_x,:)=M_full(id_homogen_x,:)+M_full(id_frame_x,:);
                M_full(:,id_homogen_x)=M_full(:,id_homogen_x)+M_full(:,id_frame_x);
                M_full(id_homogen_y,:)=M_full(id_homogen_y,:)+M_full(id_frame_y,:);
                M_full(:,id_homogen_y)=M_full(:,id_homogen_y)+M_full(:,id_frame_y);
                M_full(id_homogen_y,:)=M_full(id_homogen_y,:)+M_full(id_fluid_y,:);
                M_full(:,id_homogen_y)=M_full(:,id_homogen_y)+M_full(:,id_fluid_y);

                delete([3*iNode-2,3*iNode-1,3*iNode]) = [id_frame_x,id_frame_y,id_fluid_y];
            end

            %reducing matrix
            for i_id1=1:length(ids1)
                if ~ismember(ids1(i_id1),coupling)
                    delete(end+1)=node_array_mp1(i_id1).getDof('FRAME_DISPLACEMENT_X').getId;
                    delete(end+1)=node_array_mp1(i_id1).getDof('FRAME_DISPLACEMENT_Y').getId;
                    delete(end+1)=node_array_mp1(i_id1).getDof('FLUID_DISPLACEMENT_X').getId;
                    delete(end+1)=node_array_mp1(i_id1).getDof('FLUID_DISPLACEMENT_Y').getId;
                end
            end

            for i_id2=1:length(ids2)
                if ~ismember(ids2(i_id2),coupling)
                    delete(end+1)=node_array_mp2(i_id2).getDof('DISPLACEMENT_X').getId;
                    delete(end+1)=node_array_mp2(i_id2).getDof('DISPLACEMENT_Y').getId;
                end
            end

            % apply fixed dofs
            [~, fixedDofs] = femModel.getDofConstraints;
                if ~ isempty(fixedDofs)
                    fixedDofIds = fixedDofs.getId();
                else
                    fixedDofIds = [];
                end
                
            M=M_full;
            M(delete,:) = [];
            M(:,delete) = [];
            
            Mred=M_full;
            Mred([delete' fixedDofIds],:) = [];
            Mred(:,[delete' fixedDofIds]) = []; 
        end
        
        function [forceVector, reducedForceVector] = applyExternalForces(femModel)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            forceVector = zeros(1,nDofs);
            [~, fixedDofs] = femModel.getDofConstraints();
            
            % get the external load on the dofs
            for itDof = 1:nDofs
                if ( ~ dofs(itDof).isFixed)
                    forceVector(itDof) = dofs(itDof).getDofLoad;
                end
            end
        end
        
        function assignResultsToDofs(femModel, resultVector)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            itResult = 1;
            nResultSteps = size(resultVector,2);
            
            for itDof = 1:nDofs
                if dofs(itDof).isFixed
                    values(1:nResultSteps) = dofs(itDof).getValue(1);
                    dofs(itDof).setValue(values);
                else
                    dofs(itDof).setValue(resultVector(itResult,:));
                    itResult = itResult + 1;
                end
            end
            
        end
end
end           