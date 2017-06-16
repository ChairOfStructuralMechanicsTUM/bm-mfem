function [ output_args ] = findDuplicates( array )
%HASDUPLICATES Checks an aray for repeated values
%   Detailed explanation goes here
if (length(array) == length(unique(array)))
    output_args = false;
else
    output_args = true;

diff
end

