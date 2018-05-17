% initialize the bm-fem tools
disp('     __                    ______________  ___')
disp('    / /_  ____ ___        / ____/ ____/  |/  /')
disp('   / __ \/ __ `__ \______/ /_  / __/ / /|_/ / ')
disp('  / /_/ / / / / / /_____/ __/ / /___/ /  / /  ')
disp(' /_.___/_/ /_/ /_/     /_/   /_____/_/  /_/   ')
disp('Initializing bm-fem')

% add all folders to the search path
addpath('core');
addpath('core/assemblers');
addpath('core/solvers');
addpath('elements');
if (exist('external_libs','dir') == 7)
    addpath(genpath('external_libs'));
end
addpath('examples');
addpath('tests');
addpath('utilities');

% switch off unnecessary warnings
warning off MATLAB:handle_graphics:exceptions:SceneNode