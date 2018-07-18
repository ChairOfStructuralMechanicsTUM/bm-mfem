classdef InverseApproach < SolutionMethod
    % InverseApproach 
    %
    % Provides the abstract class for the one-, two- and three-dimensional
    % inverse approach
    %
    % Properties:
    %    numberOfBands - number of bands to compute
    %    Rref          - template of reduced matrix
    % Methods:
    %    obj = InverseApproach(model,inputInformation) - constructor
    %    R = imposedBC(obj,R,amount,index,nx,ny,nz)
    %
   properties
       numberOfBands 
       Rref
   end
   methods (Abstract=true)
       [Mred,Kred] = getReducedMatrices(obj,n); 
       f = solution(obj)
   end
   methods
       function obj = InverseApproach(model,inputInformation)
           % Description: Constructor for InverseApproach objects 
           %
           % Parameters : model -- Model object
           %              inputInformation -- Array of strings
           %
           % Return     : obj -- InverseApproach object
           obj = obj@SolutionMethod(model,inputInformation);
           obj.numberOfBands = inputInformation.nBands;
           R =  eye(length(model.M));
           sn = obj.model.slaveNodes;
           for i = obj.model.numMasterNodes:-1:1
               aRange = sn(i).a;
               bRange = sn(i).b;
               R(:,aRange :bRange) = [];
           end
           obj.Rref = R;
       end
 
       function  R = imposedBC(obj,R,miu,index,nx,ny,nz)
           % Description: Impose boundary condition on reduced matrix
           % Parameters: R     -- Reduced matrix
           %             miu   -- imposed boundary condition
           %             index -- index of node where the boundary 
           %                        condition has to be imposed on
           %             nx    -- x coordinate of the master node
           %             ny    -- y coordinate of the master node
           %             nz    -- z coordinate of the master node
           % Return:     R -- Modified reduction amtrix
           nNodes = obj.model.numNodes;
           nodes = obj.model.nodes;           
           
           for j = 1 : nNodes
               if nodes(j).coord.x==nx && ...
                  nodes(j).coord.y==ny && ...
                  nodes(j).coord.z==nz
                  rootNode = j;
                  break;
               end
           end

           for j = 1: size(R,2)
               if obj.Rref(nodes(rootNode).b,j) == 1
                   rcol = j;
                   break;
               end
           end
           
           for j = 0 : nodes(index).dof - 1
               R(nodes(index).b -j,rcol-j) = miu;
           end
       end
   end
end
