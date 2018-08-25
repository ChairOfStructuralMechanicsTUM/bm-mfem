classdef substructureFETI_DP < handle
    
   %% properties
   properties (Access=private)
       nodematrix
       
   end
   
   
   
   %% test and constructor
   methods
       
       function substructuring= substructureFETI_DP(femModel)
             if (nargin > 0)
                
            else
                error('input model is missing');
            end
       end
       
   end
   
   
  %% methods
   methods (Static)
       
       %% nodematrix
       function [nodematrix]=setupNodeMatrix(femModel,dim)
       nodearray=femModel.getAllNodes;
       nodeIdArray=zeros(dim);
       for itnod = 1:length(nodearray)
           nodeIdArray(itnod)=nodearray(itnod).getId;
       end
       
       nodematrix=zeros(dim);
       k=1;
       for i=1:dim(2)
           for j=1:dim(1)
           nodematrix(j,i)=nodeIdArray(k);
           k=k+1;
           end
       end
       end
       
       
       %% nodematrix substructuring
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
   end
end