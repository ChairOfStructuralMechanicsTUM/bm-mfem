function runTests( tests )
%RUNTESTS Summary of this function goes here
%   Detailed explanation goes here
switch tests
    case 'small'
        run(ElementTests);
        run(AssemblerTests);
        
    case 'validation'
        run(ValidationTests);
        
    case 'all'
        run(ElementTests);
        run(AssemblerTests);
        run(ValidationTests);
        
    otherwise
        disp('possible input parameters are SMALL, VALIDATION, or ALL')

end

