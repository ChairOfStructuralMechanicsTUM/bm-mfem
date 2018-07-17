function b = beam_test_nmalm_blochInv(m,n)
%2DElement with m x n nodes, 
%distance in x direction = 0.5
%distance in x direction = 1
    
    model = FemModel();
   
  for a=1:m  
    for i=1:n   
        model.addNewNode((a-1)*n+i,(i-1)*0.5,a-1); 
    end
  end %Erstellen der Knotenmatrix, mit ID, x-Koord., und y-Koord.
   
  model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

    
  for a=1:(m-1)  
    for i=1:(n-1)
        x=(a-1)*(n-1)+i; %ID des Elements
        y=(a-1)*n+i; %ID des oberen linken Knoten des Elements
        model.addNewElement('QuadrilateralElement2d4n',x,[y y+1 y+n+1 y+n]);
    end %Elemente mit ID i, bestehend aus 4 Knoten (obere und mittlere Reihe);innere Knoten jew Doppelt
  end
    
%   for a=1:m
%     model.getNode(1+(a-1)*n).fixDof('DISPLACEMENT_X');
%     model.getNode(1+(a-1)*n).fixDof('DISPLACEMENT_Y');
%   end
  
    model.getAllElements.setPropertyValue('YOUNGS_MODULUS',210000000000);
    model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
    model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
    model.getAllElements.setPropertyValue('DENSITY',7860);

%     addPointLoad(model.getNode(m*n),10,[0 -1]);
  
    obj = BlochInverse1D(model);
    assembling = SimpleAssembler(model);
    
    stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

    massMatrix = assembling.assembleGlobalMassMatrix(model);
    initialize(obj);


%     solver = SimpleSolvingStrategy(model);
%     x = solver.solve();
% 
%     step = 1;
% 
%     VerschiebungDofs = model.getDofArray.getValue(step);
% 
%     nodalForces = solver.getNodalForces(step);
% 
%     v = Visualization(model);
%     v.setScaling(10000000);
%     v.plotUndeformed
%     v.plotDeformed
    
end
