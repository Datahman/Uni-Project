function [] = PlotIXPolar(NumThetasInPi2, LinearLength, LatticeSpacing,polardata, plottype)
%  NumThetasInPi2 = 20;
%  LinearLength = 6;
% load(matlab,'-mat')
% polardata = load('matlab.mat')
% polardata
ThetaXt = zeros(1,(4 * NumThetasInPi2) + 2);
for ttt=1:(4 * NumThetasInPi2) + 2 
       ThetaXt(ttt) = ((ttt-1) * 2 * pi)/(4 * NumThetasInPi2);
end
XXX = zeros(LinearLength,(4 * NumThetasInPi2) + 2);
YYY = zeros(LinearLength,(4 * NumThetasInPi2) + 2);
CCC = zeros(LinearLength,(4 * NumThetasInPi2) + 2);
rrr = (0:LinearLength-1)' .* LatticeSpacing;  %2.25 * (0:60)'/60;
XXX = rrr*cos(ThetaXt);
YYY = rrr*sin(ThetaXt);


switch plottype 
    case 1
% % %I to I, I think.
CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =     (((polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,1))));
title('Particle density from an unpolarized source')
    case 2
% % Sx to Sz, I think.
 CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,4,2)+polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,4,1)) ./ (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,2)+polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,1)) );
title('Sx to sz spin Polarization Degree')
    case 3
% Sx to Sy, I think.
  CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,3,2)+polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,3,1)) ./ (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,2)+polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,1)) );
 title('Sx to sy spin Polarization Degree')
    case 4
% %I to Sz, I think.
 CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,4,1))) ./ real( (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,1)) );
  title('I to sz spin Polarization Degree')
    case 5
% % %I to Sy, I think.
 CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,3,1)) ./ (polardata(1:LinearLength,1:(4 * NumThetasInPi2) + 1,1,1)) );
title('I to sy spin Polarization Degree')    
    otherwise
        msgbox 'bad plot type'
end
        


% tempCCC = CCC(:,1:NumThetasInPi2);
% CCC(:,1:NumThetasInPi2) = CCC(:,(2*NumThetasInPi2)+1:3*NumThetasInPi2);
% CCC(:,(2*NumThetasInPi2)+1:3*NumThetasInPi2) = tempCCC;
% CCC(:,(4*NumThetasInPi2)+2) = CCC(:,2);
% CCC(:,(4 * NumThetasInPi2) + 1)=CCC(:,1);

CCC(:,(4 * NumThetasInPi2) + 2)=CCC(:,2);
contour(XXX,YYY,CCC,100);
colorbar;

% axis equal
% axis([-5, 5,-5,5]);
% caxis([-0.05,0.05])
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','$${\mu m}$$','Interpreter','latex','fontsize',14)
% set(gca,'XTick',[-5 0 5],'fontsize',14)
% set(gca,'YTick',[-5 0 5],'fontsize',14)
% set(gca,'YTickLabel',[])
% set(gca,'XMinorTick','off')
% set(gca,'YMinorTick','off')
% set(gca,'ticklength',5*get(gca,'ticklength'))
% set(gca,'units','centimeters')
% set(gca,'Position',[8.5 2 3.8 3.8]);

end
