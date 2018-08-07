clear;

node01 = Node(1,0,0);
node02 = Node(2,1,0);
node03 = Node(3,2,0);
node04 = Node(4,2,1);
node05 = Node(5,1,1);
node06 = Node(6,0,1);

nodeArray_homogen = [node01 node02 node05 node06];
nodeArray_porous = [node02 node03 node04 node05];

nodeArray_homogen.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
nodeArray_porous.addDof({'FRAME_DISPLACEMENT_X', 'FRAME_DISPLACEMENT_Y',...
                'FLUID_DISPLACEMENT_X', 'FLUID_DISPLACEMENT_Y'});

element_homogen = Quadrilateral2d4n(1,nodeArray_homogen);
element_porous = Porous2d4n(1,nodeArray_porous);

%Values
element_porous.setPropertyValue('FREQUENCY',100); %angenommen 100 Hz
element_porous.setPropertyValue('NUMBER_GAUSS_POINT',2);
%Frame
element_porous.setPropertyValue('FRAME_LAME_PARAMETER_LAMBDA',114400);
element_porous.setPropertyValue('FRAME_LAME_PARAMETER_MU',28600);
%Fluid
element_porous.setPropertyValue('HEAT_CAPACITY_RATIO',1.4);
element_porous.setPropertyValue('AMBIENT_FLUID_STANDARD_PRESSURE',101300); %angenommen 1013 hPa
element_porous.setPropertyValue('AMBIENT_FLUID_VISCOSITY',1.84*10^-5);
element_porous.setPropertyValue('AMBIENT_FLUID_DENSITY',1.21);
element_porous.setPropertyValue('PRANDTL_NUMBER',0.71);
%Porous
element_porous.setPropertyValue('THERMAL_CHAR_LENGTH',226*10^-6);
element_porous.setPropertyValue('POROSITY',0.9);
%Mass Matrix
element_porous.setPropertyValue('FRAME_DENSITY',300); % berechnet aus Porosität und Angabe
element_porous.setPropertyValue('TORTUOSITY',7.8);
element_porous.setPropertyValue('VISCOUS_CHAR_LENGTH',226*10^-6);
element_porous.setPropertyValue('STATIC_FLOW_RESISTIVITY',2500);
          