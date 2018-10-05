function [ duplicates_array ] = findDuplicates( array )
%HASDUPLICATES Returns repeated values in an array
%   Detailed explanation goes here
sorted_array = sort(array);
duplicates_indices = find(diff(sorted_array) == 0);
duplicates_array = sorted_array(duplicates_indices);

end

