    clear
%     clc

    io=MdpaInput('Beispiel_rechteckigerEinschluss_nachrechtsversetzt.mdpa'); %specify input file
    model = io.readModel(); %read the model

    model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
    
    a=model.getAllModelParts;
    
%     leftNodes = model.getModelPart('GENERIC_leftNodes').getNodes();
%     rightNodes = model.getModelPart('GENERIC_rightNodes').getNodes();
%     outerMaterial = model.getModelPart('GENERIC_outerMaterial').getElements();
%     innerMaterial = model.getModelPart('GENERIC_innerMaterial').getElements();
%     innerMaterialNodes = model.getModelPart('GENERIC_innerMaterial').getNodes();
%     outerMaterialNodes =  model.getModelPart('GENERIC_outerMaterial').getNodes();

    b = model.getAllNodes();
    c = model.getAllElements();
    v=Visualization(model); %set up visualization
    v.plotUndeformed()  %visualize
    
    %aluminum
    model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10); %7e10
    model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
    model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
    model.getAllElements.setPropertyValue('DENSITY',2699); %2699
    
    %structural steel
%     model.getModelPart('GENERIC_innerMaterial').getElements.setPropertyValue('YOUNGS_MODULUS',21e10); %%21e10
%     model.getModelPart('GENERIC_innerMaterial').getElements.setPropertyValue('POISSON_RATIO',0.34);
%     model.getModelPart('GENERIC_innerMaterial').getElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
%     model.getModelPart('GENERIC_innerMaterial').getElements.setPropertyValue('DENSITY',7856);   

%     steelDensity = innerMaterial(1).getPropertyValue('DENSITY');

    solver = BlochInverse1D(model);
    assembling = SimpleAssembler(model);

    stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

    massMatrix = assembling.assembleGlobalMassMatrix(model);


    initialize(solver)


    [Ksorted,Msorted] = sortKandM(solver,stiffnessMatrix,massMatrix);

    numberOfPhases = 20;

    [Kred,Mred] = reducedStiffnesAndMass (solver,Ksorted,Msorted,numberOfPhases);  

    omega = cell(numberOfPhases,1);

    nob = 10;
    [kx,miu] = propConst(solver,numberOfPhases);

    for j = 1:nob
    
        for i = 1:numberOfPhases
            omega{i,1} = solver.calcOmega(Kred{i,1},Mred{i,1},nob);
            f(j,i) = omega{i,1}(j,1)/(2*pi);

        end

        figure(2)
        plot(kx,f(j,:),'r')
        title('Dispersion curves')
        xlabel('Phase k')
        ylabel('frequenzy f')
        xlim([0 pi])
        legend(['bandnumbers: ' num2str(j)],'Location','EastOutside')
        hold on
    end
