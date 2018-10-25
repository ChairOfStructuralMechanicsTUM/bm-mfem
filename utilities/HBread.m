function [ A ] = HBread( file )
%HBREAD Reads a matrix stored in the Harwell-Boeing format to Matlab
%   Only real, symmetric, assembled (RSA) matrices are supported
fid=fopen(file) ;
if fid < 0, error('Cannot open file'); end

fgetl(fid);

header = strsplit(fgetl(fid));
PTRCRD = str2double(header{3});     %data pointers
INDCRD = str2double(header{4});     %row indices
VALCRD = str2double(header{5});     %number of data points

header = strsplit(fgetl(fid));

if ~strcmp(header{1},'RSA')
    error('Unknown matrix format. Please use hb_to_msm.m')
end

ncol = str2double(header{2});
nrow = str2double(header{3});

fgetl(fid); fgetl(fid);

COLPTR = cell2mat(textscan(fid,'%.0f',PTRCRD));
ROWIND = cell2mat(textscan(fid,'%.0f',INDCRD));
VALUES = cell2mat(textscan(fid,'%f',VALCRD));

col_ind = zeros(VALCRD,1);
a = 1;
for ii=1:length(COLPTR)-1
    n = COLPTR(ii+1) - COLPTR(ii);
    col_ind(a:a+n) = ii;
    a = a+n;
end

col_ind(end) = [];
row_ind = ROWIND;

A = sparse(row_ind, col_ind, VALUES, nrow, ncol);

fclose(fid);

end

