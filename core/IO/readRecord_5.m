%
% function nodeNum = readRecord_5()
%
%
% Function   : readRecord_5
%
% Description: This function get the nodal equivalence between the 
%              original number of node given in ANSYS and the distribution
%              within the stiffness matrix through reading record5
%
% Parameters : 
%
% Return     : nodeNum                   - array
%
function nodeNum = readRecord_5()
fid = fopen('DataAnsys/record_5.txt', 'r') ;            
if fid < 0, error('Cannot open file'); end    
for i = 1 : 10
    fgetl(fid) ;                              
end
buffer = fread(fid, Inf) ;                    
fclose(fid);
fid = fopen('DataAnsys/record_5_modified.txt', 'w')  ;  
fwrite(fid, buffer) ;                         
fclose(fid) ;
filename = 'DataAnsys/record_5_modified.txt';
delimiterIn = ' ';
A = importdata(filename,delimiterIn);   
%wenn nur ein Element in FEmodel dann ist A kein struct
if ~isstruct(A)
    B.textdata=strtrim(cellstr(num2str(A'))');
    A=B;
else
A.textdata(:,6)=[];  
end
nodeNum = [];
for i = 1 : size(A.textdata,1)               
    %for j = 1 : 5
    for j = 1 : size(A.textdata,2)      
        if isnan(str2double(A.textdata(i,j)))
            break;
        end
        nodeNum = [nodeNum str2double(A.textdata(i,j))];
    end
end
