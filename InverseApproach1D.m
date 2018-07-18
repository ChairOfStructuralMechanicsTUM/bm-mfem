classdef InverseApproach1D < InverseApproach
    % InverseApproach1D  
    %
    % Implementation of the inverse approach for one-dimensional infinite
    % structures.
    %
    % Properties:
    %    kx - wave number
    %
    % Methods:
    %    InverseApproach1D(model,inputInformation) - constructor
    %    solution(obj) - implementation of the inverse method
    %    getReducedMatrices(obj,index) - get the reduced matrices 
    %
   properties
       kx
   end
   methods
       function obj = InverseApproach1D(model,inputInformation)
           % Description: Constructor for InverseApproach1D objects 
           %
           % Parameters : model -- Model object
           %              inputInformation -- Array of strings
           %
           % Return     : obj -- InverseApproach1D object
           obj = obj@InverseApproach(model,inputInformation);
           obj.kx = linspace(1e-6,pi,obj.nSteps);
       end
       
       function f = solution(obj)
           % Description: Solve the one-dimensional model through the 
           %              inverse approach
           %
           % Parameters : obj
           %
           % Return     : f  - frequency values
           %             phi - modes
           nSteps = obj.nSteps;
           nb = obj.numberOfBands;
           f =  zeros(nb,nSteps);
           if obj.parallelProcessing == 1
               parfor j = 1 : nSteps
                   [Mred, Kred] = obj.getReducedMatrices(j);
                   omega2 = eigs(Kred,Mred,nb,'sm');
                   f(:,j)= sqrt(abs(omega2))/(2*pi);
               end
           else
               for j = 1 : nSteps
                   %j
                   [Mred,Kred] = obj.getReducedMatrices(j);
                   omega2 = eigs(Kred,Mred,nb,'sm');
                   f(:,j)= sqrt(abs(omega2))/(2*pi);
               end
           end
       end
       
       function [Mred,Kred] = getReducedMatrices(obj,index)
           % Description: Function to compute the mass and stiffness
           %              matrices according to certain wave number
           %
           % Parameters : obj
           %              index
           %
           % Return     : Mred - Reduced mass matrix
           %              Kred - Reduced stiffness matrix
           R = obj.Rref;
           miu = exp(1i*obj.kx(index));
           sNodes = obj.model.slaveNodes;
           for i = 1 : obj.model.numMasterNodes
               y = sNodes(i).coord.y;
               dof = obj.model.nodes(i).dof;
               if dof ~= 1
                   ind = sNodes(i).ind;              
                   R = obj.imposedBC(R,miu,ind,0,y,0);
               end
           end
           Mred = R'*(obj.model.M*R);
           Kred = R'*(obj.model.K*R);
       end
   end
end
