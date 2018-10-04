function runTests( tests )
%RUNTESTS Summary of this function goes here
%   Detailed explanation goes here
if nargin ~= 1
    tests = 'all';
end

switch tests
    case 'small'
        run(AssemblerTests);
        run(ElementTests);
        run(SolverTests);
        run(VariousTests);
        
    case 'validation'
        run(ValidationTests);
        
    case 'all'
        run(ElementTests);
        run(AssemblerTests);
        run(SolverTests);
        run(ValidationTests);
        run(VariousTests);
        
    otherwise
        disp('possible input parameters are SMALL, VALIDATION, or ALL')

end

