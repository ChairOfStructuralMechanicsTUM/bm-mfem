%
% function A = readRestrictions()
%
% Function   : readRestrictions
%
% Description: This function get the restrictions on the coordinates from
%              ANSYS
%
% Parameters :
%
% Return     : A array with restrictions
function A = readRestrictions()
fid=fopen('DataAnsys/nodeRest.txt') ;                   % the original file
fidd=fopen('DataAnsys/nodeRest_modified.dat','w') ;     % the new file
if fid < 0, error('Cannot open file'); end    % Check for an error
for j = 1 : 13
    fgetl(fid) ;                              % Read/discard line.
end
while ~feof(fid)  % reads the original till last line
  tline=fgets(fid);  
     if isspace(tline) 
         for j = 1 : 9
             fgetl(fid) ;                     % Read/discard line.
         end
     else
       fwrite(fidd,tline) ;
     end
end
fclose all ;
filename = 'DataAnsys/nodeRest_modified.dat';
delimiterIn = ' ';
A = importdata(filename,delimiterIn);        % Get data in matlab

if ~isempty(A)
    A.textdata(:,2) = [];
    nodeRest = str2double(A.textdata);
    A = sort(unique(nodeRest),'ascend');
else
    A = [];
end

