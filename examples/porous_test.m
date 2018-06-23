clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,2,2);
node04 = Node(4,0,2);

nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'FRAME_DISPLACEMENT_X', 'FRAME_DISPLACEMENT_Y'
                'FLUID_DISPLACEMENT_X', 'FLUID_DISPLACEMENT_Y'});

element = Porous2d4n(1,nodeArray);

%Values
element.setPropertyValue('FREQUENCY',100); %angenommen 100 Hz
element.setPropertyValue('NUMBER_GAUSS_POINT',2);
%Frame
element.setPropertyValue('FRAME_LAME_PARAMETER_LAMBDA',1144);
element.setPropertyValue('FRAME_LAME_PARAMETER_MU',286);
%Fluid
element.setPropertyValue('HEAT_CAPACITY_RATIO',1.4);
element.setPropertyValue('AMBIENT_FLUID_PRESSURE',1013); %angenommen 1013 hPa
element.setPropertyValue('AMBIENT_FLUID_VISCOSITY',1.84*10^-5);
element.setPropertyValue('AMBIENT_FLUID_DENSITY',1.21);
element.setPropertyValue('PRANDTL_NUMBER',0.71);
%Porous
element.setPropertyValue('THERMAL_CHAR_LENGTH',226);
element.setPropertyValue('POROSITY',0.9);
%Mass Matrix
element.setPropertyValue('FRAME_DENSITY',300); % berechnet aus Porosität und Angabe
element.setPropertyValue('TORTUOSITY',7.8);
element.setPropertyValue('VISCOUS_CHAR_LENGTH',226);
element.setPropertyValue('STATIC_FLOW_RESISTIVITY',25);

[K_ss,K_sf,K_fs,K_ff] = element.computeLocalStiffnessMatrix;
[M_ss,M_sf,M_fs,M_ff] = element.computeLocalMassMatrix;

