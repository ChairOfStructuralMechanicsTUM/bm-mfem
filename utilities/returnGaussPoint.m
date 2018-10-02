function [ w,g ] = returnGaussPoint( elementType, number )
%returnGaussPoint Returns evaluation point g and weighting w depending on
%number of GaussPoints

% IMPORTANT REMARK: Please add the input parameter ElementType in your
% function returnGaussPoint, e.g. Quadrilateral Element.

% if-clause required as different descriptions of Gauss-Points for
% different elements are required.


%% QUADRILATERAL and HEXAHEDRAL elements
if (strcmp(elementType,'HexahedronElement3d8n') | strcmp(elementType,'QuadrilateralElement') | strcmp(elementType,'QuadrilateralElement2d4n') ...
        | strcmp(elementType,'DiscreteKirchhoffElement3d4n') | strcmp(elementType,'ReissnerMindlinElement3d4n') | strcmp(elementType,'ShellElement3d4n') |strcmp(elementType,'~') ) 
    if number == 0
        msg = 'ReturnGaussPoint: Please specify a number > 0.';
        e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
        throw(e);
    elseif number == 1
        g=0;
        w=2;
    elseif number == 2
        g=[(-1/sqrt(3)), (1/sqrt(3))];
        w=[1,1];
    elseif number == 3
        g=[-sqrt(3/5), 0, sqrt(3/5)];
        w=[5/9, 8/9, 5/9];
    elseif number == 4
        g=[-sqrt((3/7)+((2/7)*sqrt(6/5))),-sqrt((3/7)-((2/7)*sqrt(6/5))),sqrt((3/7)-((2/7)*sqrt(6/5))),sqrt((3/7)+((2/7)*sqrt(6/5)))];
        w=[(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];
    elseif number == 5
        g=[-1/3*sqrt(5+(2*sqrt(10/7))), -1/3*sqrt(5-(2*sqrt(10/7))), 0, 1/3*sqrt(5-(2*sqrt(10/7))), 1/3*sqrt(5+(2*sqrt(10/7)))];
        w=[(322-(13*sqrt(70)))/900, (322+(13*sqrt(70)))/900, 128/225, (322+(13*sqrt(70)))/900, (322-(13*sqrt(70)))/900];
    else
        msg = 'ReturnGaussPoint: Please specify a number between 1 and 5.';
        e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
        throw(e);
    end

    
%% TETRAHEDRAL elements
% to proof: useable for hexahedral elements???       
elseif (strcmp(elementType,'TetrahedronElement3d4n'))
    if number == 0
        msg = 'ReturnGaussPoint: Please specify a number > 0.';
        e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
        throw(e);
    elseif number == 1
        g=[1/4];
        w=[1];
    elseif number == 2
        msg = 'ReturnGaussPoint: Please specify a number equal to 1 or 4.';
        e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
        throw(e); 
    elseif number == 2
        msg = 'ReturnGaussPoint: Please specify a number equal to 1 or 4.';
        e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
        throw(e);
    elseif number == 4
        g=[(5-sqrt(5))/20, (5-sqrt(5))/20, (5-sqrt(5))/20, (5+3*sqrt(5))/20];
        w=[1/4, 1/4, 1/4, 1/4];
    else
        msg = 'ReturnGaussPoint: Please specify a number equal to 1 or 4.';
        e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
        throw(e);
    end
    
    
end
    
end
