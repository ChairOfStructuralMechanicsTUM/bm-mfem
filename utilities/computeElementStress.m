function [ stressValue, prStressDir, element_connect ] = computeElementStress( elementArray, nodeArray, step )
%COMPUTEELEMENTSTRESS stress in an array of elements
%   [STRESS, PRSTRESSDIR, EC] = computeElementStress(ELEMENTS,NODES)
%   returns element connectivity EC, the matrix STRESS containing:
%       STRESS(1,:) = stress in xx
%       STRESS(2,:) = stress in yy
%       STRESS(3,:) = stress in xy
%       STRESS(4,:) = first principal stress
%       STRESS(5,:) = second principal stress
%       STRESS(6,:) = von-Mises stress
%   and the matrix PRSTRESSDIR containing:
%       PRSTRESSDIR(1,:,:) = direction of principal stress I
%       PRSTRESSDIR(2,:,:) = direction of principal stress II
%
%   [STRESS, PRSTRESSDIR, EC] = computeElementStress(ELEMENTS,NODES,STEP)
%   returns STRESS and PRSTRESSDIR for the selected STEP.
%
%   This function can only be used if all ELEMENTS are of the same type.

if nargin == 2
    step = 1;
end

nElements = length(elementArray);
nNodes = length(nodeArray);
nElementalNodes = length(elementArray(1).getNodes());

element_connect = zeros(nElements,nElementalNodes);
stressValue = zeros(6,nNodes);
sigma_xx = zeros(nElements,nElementalNodes);
sigma_yy = zeros(nElements,nElementalNodes);
sigma_xy = zeros(nElements,nElementalNodes);
smooth_sigma_xx = zeros(1,nNodes);
smooth_sigma_yy = zeros(1,nNodes);
smooth_sigma_xy = zeros(1,nNodes);
prin_I = zeros(1,nNodes);
prin_II = zeros(1,nNodes);
vm_stress = zeros(1,nNodes);
prStressDir = zeros(2,nNodes,2);

for i = 1:nElements
    element_connect(i,:) = elementArray(i).getNodes.getId();
    stressPoints = elementArray(i).stressPoints();
    EModul = elementArray(i).getPropertyValue('YOUNGS_MODULUS');
    prxy = elementArray(i).getPropertyValue('POISSON_RATIO');
    % Moment-Curvature Equations
    D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
    % Material Matrix D
    D = D * EModul / (1-prxy^2);
    
    for j = 1:nElementalNodes
        if size(stressPoints,2) == 2 %quadrilateral elements
            [~, ~, B, ~] = elementArray(i).computeShapeFunction(stressPoints(j,1),stressPoints(j,2));
        elseif size(stressPoints,2) == 3 %triangular elements
            [~, ~, B, ~] = elementArray(i).computeShapeFunction(stressPoints(j,:));
        else
            msg = ['ComputeElementStress: Invalid number of stress evaluation', ...
                ' points.'];
            e = MException('MATLAB:bm_mfem:stressEvaluation',msg);
            throw(e);
        end
        displacement_e = getValuesVector(elementArray(i),step);
        displacement_e = displacement_e';
        strain_e = B * displacement_e;
        stress_e = D * strain_e;
        
        % elementwise stress calculation
        sigma_xx(i,j) = stress_e(1);
        sigma_yy(i,j) = stress_e(2);
        sigma_xy(i,j) = stress_e(3);
        
    end
end

for k = 1 : nNodes
    [I,J] = find(element_connect == k);
    
    sum_sigma_xx = 0;
    sum_sigma_yy = 0;
    sum_sigma_xy = 0;
    
    for l = 1: length(I)
        sum_sigma_xx = sum_sigma_xx + sigma_xx(I(l),J(l));
        sum_sigma_yy = sum_sigma_yy + sigma_yy(I(l),J(l));
        sum_sigma_xy = sum_sigma_xy + sigma_xy(I(l),J(l));
    end
    
    smooth_sigma_xx(k) = sum_sigma_xx/length(I);
    smooth_sigma_yy(k) = sum_sigma_yy/length(I);
    smooth_sigma_xy(k) = sum_sigma_xy/length(I);
    
    stress_ele = [smooth_sigma_xx(k) smooth_sigma_xy(k);
        smooth_sigma_xy(k) smooth_sigma_yy(k)];
    
    [vec, lambda] = eig(stress_ele);
    prStressDir(1,k,:) = vec(:,1);
    prStressDir(2,k,:) = vec(:,2);
    
    prin_I(k) = lambda(1,1);
    prin_II(k) = lambda(2,2);
    
    vm_stress(k) = sqrt(prin_I(k).^2 + prin_II(k).^2 - prin_I(k) * prin_II(k));
end

stressValue(1,:) = smooth_sigma_xx;
stressValue(2,:) = smooth_sigma_yy;
stressValue(3,:) = smooth_sigma_xy;
stressValue(4,:) = prin_I;
stressValue(5,:) = prin_II;
stressValue(6,:) = vm_stress;

end

