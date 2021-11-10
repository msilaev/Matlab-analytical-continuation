function ggA = OPfun(Om1,D,h,b)

%ggA= (Om1.*Om2 - D^2)./(s1.*s2);

ggA=(i/4)*(sqrt(-1)*(-4)).*(D+D.*(b.^2+D.^2+(-1).*h.^2+Om1.^2).*((-1).*b.^2+ ...
  D.^2+(-1).*h.^2+Om1.^2+(sqrt(-1)*(-2)).*(h.^2.*Om1.^2+b.^2.*(D.^2+ ...
  Om1.^2)).^(1/2)).^(-1/2).*((-1).*b.^2+D.^2+(-1).*h.^2+Om1.^2+( ...
  sqrt(-1)*2).*(h.^2.*Om1.^2+b.^2.*(D.^2+Om1.^2)).^(1/2)).^(-1/2)).* ...
  (((-1).*b.^2+D.^2+(-1).*h.^2+Om1.^2+(sqrt(-1)*(-2)).*(h.^2.* ...
  Om1.^2+b.^2.*(D.^2+Om1.^2)).^(1/2)).^(1/2)+((-1).*b.^2+D.^2+(-1).* ...
  h.^2+Om1.^2+(sqrt(-1)*2).*(h.^2.*Om1.^2+b.^2.*(D.^2+Om1.^2)).^( ...
  1/2)).^(1/2)).^(-1);

  % s1=sqrt(  D^2 + Om1.^2 ) ;
% s2=sqrt(  D^2 + Om2.^2 ) ;
% ggA= ( Om1.*Om2 + D^2)./(s1.*s2);

end

       







% =========================================================================


