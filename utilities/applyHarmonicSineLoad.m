function applyHarmonicSineLoad( nodeArray, load, direction, excitationFrequency, time )
%APPLYHARMONICSINELOAD Applies a harmonic sine load
%   nodeArray: nodes, where the load is imposed
%   load: the magnitude of the applied load
%   direction: the direction in [x,y,z]
%   excitationFrequency: the excitation frequency of the load in rad/s
%   time: the current time

%modulus wrt time:
timeLoad = load * sin(excitationFrequency * time);
addPointLoad(nodeArray, timeLoad, direction);

end

