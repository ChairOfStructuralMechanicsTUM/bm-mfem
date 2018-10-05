function [ w,g ] = returnGaussPoint( number )
%returnGaussPoint Returns evaluation point g and weighting w depending on
%number of GaussPoints for quadrilateral and hexahedral elements

if number == 1
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

end


    
    

    

