function [ w,g ] = returnGaussPointTetrahedron( number )
%returnGaussPointTetrahedron Returns evaluation point g and weighting w
%for tetrahedral elements depending on number of GaussPoints

if number == 0
    msg = 'ReturnGaussPoint: Please specify a number > 0.';
    e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
    throw(e);
elseif number == 1
    g=1/4;
    w=1;
elseif number == 4
    g=[(5-sqrt(5))/20, (5-sqrt(5))/20, (5-sqrt(5))/20, (5+3*sqrt(5))/20];
    w=[1/4, 1/4, 1/4, 1/4];
else
    msg = 'ReturnGaussPoint: Please specify a number equal to 1 or 4.';
    e = MException('MATLAB:bm_mfem:invalidNumberOfGaussPoints',msg);
    throw(e);
end

end