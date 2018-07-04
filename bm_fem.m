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
addpath(genpath('examples'));
addpath('tests');
addpath(genpath('utilities'));

% check, if all neccessary toolboxes are installed
requiredToolboxes = {'MATLAB','Mapping Toolbox'};
v = ver;
[installedToolboxes{1:length(v)}] = deal(v.Name);
missing = setdiff(requiredToolboxes,installedToolboxes);
if ~isempty(missing)
    error('required %s is not available\n',missing{:})
end
clear v installedToolboxes requiredToolboxes missing

% switch off unnecessary warnings
warning off MATLAB:handle_graphics:exceptions:SceneNode