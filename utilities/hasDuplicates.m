function [ output_args ] = hasDuplicates( array )
%HASDUPLICATES Checks an aray for repeated values
%   Detailed explanation goes here
if (length(array) == length(unique(array)))
    output_args = false;
else
    output_args = true;


end

