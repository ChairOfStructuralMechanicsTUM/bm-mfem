function [ w,g ] = returnGaussPoint( number )
%returnGaussPoint Returns evaluation point g and weighting w depending on
%number of GaussPoints
%   Detailed explanation goes here
  if number == 1
      g=0;
      w=2;
  elseif number == 2
      g=[];
      w=[1,1];
  elseif number == 3
      g=[];
      w=[5/9, 8/9, 5/9];
  elseif number == 4
      g=[];
      w=[(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];
  elseif number == 5
      g=[];
      w=[(322-(13*sqrt(70)))/900, (322+(13*sqrt(70)))/900, 128/225, (322+(13*sqrt(70)))/900, (322-(13*sqrt(70)))/900];
  end
  
end
