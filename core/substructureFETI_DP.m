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
       function [nodematrix]=setupNodematrix(femModel)
       nodeArray = getAllNodes(femModel)
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
   end
end