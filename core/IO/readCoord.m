% function A = readCoord()
%
% Function   : readCoord
%
% Description: This function gets the nodes coordinates from txt files
%              created in ANSYS
%
% Parameters : 
%
% Return     : A                   - matrix with nodal information
%
function A = readCoord()
fid=fopen('DataAnsys/nodeCoor.txt') ;                   
fidd=fopen('DataAnsys/nodeCoor_modified.dat','w') ;     
if fid < 0, error('Cannot open file'); end 
% Discard some line to read the data from the txt files
for j = 1 : 13
    fgetl(fid) ;                              
end

while ~feof(fid) 
  tline=fgets(fid);  
     if isspace(tline) 
         for j = 1 : 9
             fgetl(fid) ;                     
         end
     else
       fwrite(fidd,tline) ;
     end
end

fclose all ;
filename = 'DataAnsys/nodeCoor_modified.dat';
delimiterIn = ' ';
% Get data in matlab
A = importdata(filename,delimiterIn); 