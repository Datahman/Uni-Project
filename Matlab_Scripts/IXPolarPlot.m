% todo: I can check the renormalization in  the holes article just by calculating the diffuson decay times numerically and comparing to the analytical formulas.  Is there a way to get the precession constants numerically, i.e. by selecting a good subspace? If I find the expected minima I can compare their positions to the analytical values, and then calculate the linear in q and quadratic in q moments of the diffuson around those minima.  These moments can be reduced to the subspace of the long-lived mode, I think, and then compared to the analytical formulas.
% 
% It might also be interesting to extend the holes article to the EY regime, using numerical evaluation
% 
% I can check the polariton article more easily because here the focus is on the diffusive approximation and its validity.  However it will be necessary to analytically expand to first order in Eso/Ef before doing the numerical integrals.  After doing this, just produce the same graphs as in the polariton article, and compare.
% 

% todo: hbarovertau must always be much smaller than Ef - otherwise the
% system will localize.  Check this.  In consequence, if Eso is similar to
% Ef then only the EY regime is allowed, not the DP regime.  Related
% question: how big can we make Ef before the exciton unbinds? I thought
% about this before.... Littlewood's article?  Szymanska and Littlewood.
% Szymanska and Littlewood PRB 67 193305 put the binding energy of the
% indirect exciton between 3 and 4 meV, or 3000-4000 micro eV.

%polar6red = IXPolarPlot(4,0,0.5, 20,60,100,8,8,100,0.5,6,4, 5);
%polar6 = IXPolarPlot(4,0,0.5, 20,60,100,8,8,100,0.5,6,10, 5); t= 744.376395 seconds.
%polar50 = IXPolarPlot(4,0,0.5, 20,60,100,8,8,100,0.5,50,40, 5);

%polar10 = IXPolarPlot(4,0,0.5, 20,60,100,8,8,100,0.5,10,10, 3);
%Elapsed time is 330.892542 seconds. 328.515990 if I comment out the expm. 
%and fft. 251.244147 after some optimization.  Then 285. Then 213 without
%the array update. 3 if all math is removed. 23 if only the vf computation.
% So: 3 overhead, 20 vf, 190 other arithmetic, 72 array update. I checked
% this again and it worked out.  Then I put everything together and it was
% 298, and then again  and it was 283.

% plots of Fermi surfaces kf as a function of theta at  B=0,1,3,5,7, Ef=20,35,50
% kf, vf, deltaE's, D's as a function of B, and as function of Ef
% Ef > 15 microeV
% 0 < B < 7 Tesla
% I believe that we have butchered the ballistic physics and can not expect reasonable behavior at small distance and length scales.
% The lattice spacing should be the same or 1/2 of the scattering length. 
% Equivalently the diffuson dispersion should be cutoff at at the scattering energy scale, but I think that is true anyway, I don't think that the correct diffuson can exceed 1.  I haven't proved this though.
% The sample size should be 5 or 10 times the spin decay length, which is like saying that the energy scale of the .
% 
% PRB 195309 sets kf = 15.3, with eigenvalues at 30, 6, -13 micro eV, and precession lengths  are shown around 2 micro meters.  The source has a radius of 4 micro m, and the sample size is 20x20 micro m  Decay takes too long, like 10 ns to 10 micro s.
% 
% The PRL 110 246403 and supporting material does spatial averaging over 1.5 microm, and the ring radius is 4 micro m. They populate the ring with T=0.1 k=0 excitons and then watch what happens. They show the linear polarization of the bright excitons. The plotted area is 30 x 30 micro m. 

function[spintheta, brightspinxy] = IXPolarPlot(Norbitals,BMagnetic,deltab,Ef, kmax,NumKGridPoints,MaxNumkf,NumValueTheta,hbarOverTau, LatticeSpacingOverUnitLength, LinearLength, NumThetasInPi2, DeltaTOverTau, UseDiffExpansion)
tic

         brightspinxy = zeros(LinearLength, LinearLength,4,4);


DinverseArray = zeros(LinearLength,LinearLength,Norbitals * Norbitals, Norbitals * Norbitals); % Diffuson operator 16x16
DOS = 0; 
kfatzero=zeros(MaxNumkf);
KfAtThetaPi2=zeros(MaxNumkf);
RhoN = zeros(Norbitals, Norbitals);

BrightIXPauliVectors = zeros(4, Norbitals * Norbitals);
% The identity matrix for bright excitons, represented as a vector.
BrightIXPauliVectors(1,1) = 1;
BrightIXPauliVectors(1,6) = 1;
% The sigma_x matrix for bright excitons, represented as a vector.
BrightIXPauliVectors(2,2) = 1;
BrightIXPauliVectors(2,5) = 1;
% The sigma_y matrix for bright excitons, represented as a vector.
% todo: there is a possible factor of minus one error in sigma_y,
%     or equivalently a confusion about the transverse of 
% todotodo: why is this not multiplied by i?
BrightIXPauliVectors(3,2) = 1;
BrightIXPauliVectors(3,5) = -1;
% The sigma_z matrix for bright excitons, represented as a vector.
BrightIXPauliVectors(4,1) = 1;
BrightIXPauliVectors(4,6) = -1;

   % %used for polar plot
% % StatsspintimeBT20xi0p25ThetaB112Cut1SO1
% % xi = 0.25, ThetaB = 112, t = 20, time
% xi = 0.25;
% ThetaB = 112.5 * pi / 180;
% tau = 4;
% timett = 20;
% spintimeBT20xi0p25ThetaB112Cut1SO1 = polariton(timett/tau, 0.2, 256, 1.0, xi, ThetaB, 0.25/pi, 200/tau, NumThetas, WhichCutoff, WhichSO);
% 
    spinxy = zeros(LinearLength,LinearLength, Norbitals * Norbitals, Norbitals * Norbitals);
    spintheta =  zeros(LinearLength, (4*NumThetasInPi2)+1,Norbitals, Norbitals);
    tmppk = zeros(Norbitals * Norbitals, Norbitals * Norbitals);
    
   for t=1:NumThetasInPi2  
       t
       ThetaX = ((t-1) * pi)/(2 * NumThetasInPi2);
kspacegrid = excitonkspace(LatticeSpacingOverUnitLength, LinearLength, ThetaX);

% this function calculates the scalar and matrix DOS, and Dinverse, after integration
% over the Fermi surface.
% It also calculates the Fermi momenta at theta = 0 and theta = pi/2. 
% It does this at particular values of Ef, tau, and the Hamiltonian parameters.  
DiffExpCoeffs = zeros(6,1,1);
if ~UseDiffExpansion
%function [IntegratedDOS, kfAtThetaZero, KfAtThetaPi2, IntegratedRhoN, IntegratedDinverse, DiffExpCoeffs] =  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumEf,NumValueTheta, hbarOverTau, kspacegrid, LinearLength)
 [DOS,kfatzero,KfAtThetaPi2,RhoN, DinverseArray,DiffExpCoeffs]=  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumkf,NumValueTheta,hbarOverTau,kspacegrid, LinearLength);
else
    kspacegriddummy = excitonkspace(1,1,0);
 [DOS,kfatzero,KfAtThetaPi2,RhoN, DinverseArray,DiffExpCoeffs]=  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumkf,NumValueTheta,hbarOverTau,kspacegriddummy, 1);
 SimplifiedDiffxx = diag(squeeze(DiffExpCoeffs(4,:,:)));
 SimplifiedDiffyy = diag(squeeze(DiffExpCoeffs(5,:,:)));
 DiffConstMax = -1;
 for i=1:Norbitals*Norbitals
     if abs(SimplifiedDiffxx(i)) > DiffConstMax
         DiffConstMax = abs(SimplifiedDiffxx(i));
     end
     if abs(SimplifiedDiffyy(i)) > DiffConstMax
         DiffConstMax = abs(SimplifiedDiffyy(i));
     end
 end
 DinverseArray = zeros(LinearLength,LinearLength,Norbitals * Norbitals, Norbitals * Norbitals); % Diffuson operator 16x16
    for qx=1:LinearLength
        for qy=1:LinearLength
            iqx = 1i * kspacegrid(qx,qy,1) *  kspacegrid(qx,qy,2);
            iqy = 1i * kspacegrid(qx,qy,1) *  kspacegrid(qx,qy,3);
            qlsqr = 2 * DiffConstMax * kspacegrid(qx,qy,1) * kspacegrid(qx,qy,1);
            Dtmp = DiffExpCoeffs(1,:,:) + (iqx * DiffExpCoeffs(2,:,:)) + (iqy * DiffExpCoeffs(3,:,:));
                + (iqx * iqx * DiffExpCoeffs(4,:,:)) + (2 * iqy * iqy * DiffExpCoeffs(5,:,:)) + (iqx * iqy * DiffExpCoeffs(6,:,:));
   
                Dtmp2 = squeeze(DiffExpCoeffs(1,:,:));
            for nn = 1:Norbitals*Norbitals
                Dtmp2(nn,nn) = Dtmp2(nn,nn) + (iqx * iqx * SimplifiedDiffxx(nn)) + (iqy * iqy * SimplifiedDiffyy(nn));
            end
            if qlsqr < 0.5
                DinverseArray(qx,qy,:,:) = squeeze(Dtmp);
            elseif qlsqr > 2.0
                DinverseArray(qx,qy,:,:) = Dtmp2;
            else
                DinverseArray(qx,qy,:,:) = (squeeze(Dtmp) * (2.0 - qlsqr)/1.5) + (Dtmp2 * (qlsqr - 0.5)/1.5);
            end
        end
    end
end

       for i=1:LinearLength
           for j=1:LinearLength
                       tmppk = squeeze(DinverseArray(i,j,:,:));
               %DinverseArray(i,j,:,:) = tmppk;
               DinverseArray(i,j,:,:) = expm(-tmppk * DeltaTOverTau);
            end
       end
       
        %spinxy = DinverseArray;
        spinxy = fft2(DinverseArray);
         brightspinxy = zeros(LinearLength, LinearLength,4,4);
       for i=1:LinearLength
           for j=1:LinearLength
                brightspinxy(i,j,:,:) = BrightIXPauliVectors * squeeze(spinxy(i,j,:,:)) * BrightIXPauliVectors';
            end
       end
      
        
      
          spintheta(:, t,:,:) = brightspinxy(1,:,:,:);
        spintheta(:, t + NumThetasInPi2,:,:) = brightspinxy(:,1,:,:);


          for i=1:LinearLength-1
            spintheta(i+1, t + (2 * NumThetasInPi2),:,:) = brightspinxy(1,LinearLength+1-i,:,:);
            spintheta(i+1, t + (3 * NumThetasInPi2),:,:) = brightspinxy(LinearLength+1-i,1,:,:);
        end
            spintheta(1, t + (2 * NumThetasInPi2),:,:) = brightspinxy(1,1,:,:);
            spintheta(1, t + (3 * NumThetasInPi2),:,:) = brightspinxy(1,1,:,:);
%            spintheta((4 * NumThetasInPi2) + 1, t,:,:) = spintheta(1, t,:,:);
        end
 
temp = spintheta(:,1:NumThetasInPi2,:,:);
spintheta(:,1:NumThetasInPi2,:,:) = spintheta(:,(2*NumThetasInPi2)+1:3*NumThetasInPi2,:,:);
spintheta(:,(2*NumThetasInPi2)+1:3*NumThetasInPi2,:,:) = temp;
spintheta(:,(4*NumThetasInPi2)+2,:,:) = spintheta(:,2,:,:);
spintheta(:,(4 * NumThetasInPi2) + 1,:,:)=spintheta(:,1,:,:);

toc
end

% % % Sujeong's 16-color palette color_vs.docx
% % 244,119,193
% % 233,37,100
% % 255,0,0
% % 193,68,18
% % 249,157,57
% % 163,161,161
% % 255,242,0
% % 130,123,0
% % 57,181,74
% % 0,116,107
% % 0,174,228
% % 0,84,166
% % 69,99,189
% % 31,58,104
% % 146,39,143
% % 0,0,0
% % 
% % %Sujeong's 8-color palette color_vs.docx
% % 244,119,193
% % 255,0,0
% % 249,157,57
% % 255,242,0
% % 57,181,74
% % 0,174,228
% % 69,99,189
% % 146,39,143
% % 
% % % Sujeong's 6-color palette color_vs.docx
% % 244,119,193
% % 255,0,0
% % 242,101,34
% % 186,221,22
% % 0,118,163
% % 0,0,0
% % 
% % % Sujeong's 4-color palette color_vs.docx
% % 255,0,0
% % 242,101,34
% % 186,221,22
% % 0,118,163
% % 
% 
%  set(0,'DefaultTextFontname', ' Times New Roman ')
%    set(0,'DefaultAxesFontName', ' Times New Roman ')
% 
% ThetaXt = zeros(1,162);
% for ttt=1:162 
%        ThetaXt(ttt) = ((ttt-1) * 2 * pi)/(160);
% end
% XXX = zeros(31,162);
% YYY = zeros(31,162);
% 
% str11 = sprintf('(a)');
% str12 = sprintf('(b)');
% str13 = sprintf('(c)');
% str14 = sprintf('(d)');
% 
% 
% 
% subplot(2,2,1);
% hold off
% rrr = (0:30)' * 0.2 /sqrt(20/(2 * 4));  %2.25 * (0:60)'/60;
% %rrr = 0.75 * 3 * (0:30)'/30;
% XXX = rrr*cos(ThetaXt);
% YYY = rrr*sin(ThetaXt);
% CCC = zeros(31,162);
% tempCCC = zeros(31,40);
% CCC(1:31,1:161) =  real( (spintimeBT20xi0p25ThetaB112Cut1SO1(1:31,:,4,2)+spintimeBT20xi0p25ThetaB112Cut1SO1(1:31,:,4,1)) ./ (spintimeBT20xi0p25ThetaB112Cut1SO1(1:31,:,1,2)+spintimeBT20xi0p25ThetaB112Cut1SO1(1:31,:,1,1)) );
% tempCCC = CCC(:,1:40);
% CCC(:,1:40) = CCC(:,81:120);
% CCC(:,81:120) = tempCCC;
% CCC(1:31,162) = CCC(1:31,2);
% %CCC =  real( (spintimeBT20xi0p25ThetaB112Cut1SO2(1:31,:,4,2)+spintimeBT20xi0p25ThetaB112Cut1SO2(1:31,:,4,1)) ./ (spintimeBT20xi0p25ThetaB112Cut1SO2(1:31,:,1,2)+spintimeBT20xi0p25ThetaB112Cut1SO2(1:31,:,1,1)) );
% contour(XXX,YYY,CCC,20);
% hold on
% shading flat
% text(-4.75,4.4,str11, 'fontsize',14) 
% colorbar;
% axis equal
% axis([-5, 5,-5,5]);
% caxis([-0.5,0.5])
% xlhand = get(gca,'xlabel');
% ylhand = get(gca,'ylabel');
% set(ylhand,'string','$${y/\sqrt{Dt}}$$','Interpreter','latex','fontsize',14)
% set(gca,'XTick',[-5 0 5],'fontsize',14)
% set(gca,'YTick',[-5 0 5],'fontsize',14)
% set(gca,'XTickLabel',[' ',' ',' '])
% set(gca,'XMinorTick','off')
% set(gca,'YMinorTick','off')
% set(gca,'ticklength',5*get(gca,'ticklength'))
% set(gca,'units','centimeters')
% set(gca,'Position',[2 7.0 3.8 3.8]);
% %set(gca,'Position',[2 6.6 3.8 3.8]);
% 
% 
% 
% subplot(2,2,2);
% hold off
% rrr = (0:60)' * 0.2 /sqrt(80/(2 * 4));  %2.25 * (0:60)'/60;
% XXX = rrr*cos(ThetaXt);
% YYY = rrr*sin(ThetaXt);
% CCC = zeros(61,162);
% tempCCC = zeros(61,40);
% CCC(1:61,1:161) =  real( (spintimeBT80xi0p25ThetaB112Cut1SO1(1:61,:,4,2)+spintimeBT80xi0p25ThetaB112Cut1SO1(1:61,:,4,1)) ./ (spintimeBT80xi0p25ThetaB112Cut1SO1(1:61,:,1,2)+spintimeBT80xi0p25ThetaB112Cut1SO1(1:61,:,1,1)) );
% tempCCC = CCC(:,1:40);
% CCC(:,1:40) = CCC(:,81:120);
% CCC(:,81:120) = tempCCC;
% CCC(1:61,162) = CCC(1:61,2);
% %CCC =  real( (spintimeBT80xi0p25ThetaB112Cut1SO2(1:81,:,4,2)+spintimeBT80xi0p25ThetaB112Cut1SO2(1:81,:,4,1)) ./ (spintimeBT80xi0p25ThetaB112Cut1SO2(1:81,:,1,2)+spintimeBT80xi0p25ThetaB112Cut1SO2(1:81,:,1,1)) );
% contour(XXX,YYY,CCC,20);
% hold on
% shading flat
% text(-4.75,4.4,str12, 'fontsize',14) 
% colorbar;
% axis equal
% axis([-5, 5,-5,5]);
% caxis([-0.05,0.05])
% xlhand = get(gca,'xlabel');
% ylhand = get(gca,'ylabel');
% set(gca,'XTick',[-5 0 5],'fontsize',14)
% set(gca,'YTick',[-5 0 5],'fontsize',14)
% set(gca,'YTickLabel',[])
% set(gca,'XTickLabel',[' ',' ',' '])
% set(gca,'XMinorTick','off')
% set(gca,'YMinorTick','off')
% set(gca,'ticklength',5*get(gca,'ticklength'))
% set(gca,'units','centimeters')
% set(gca,'Position',[8.5 7.0 3.8 3.8]);
% %set(gca,'Position',[8.2 6.6 3.8 3.8]);
% 
% 
% subplot(2,2,3);
% hold off
% rrr = (0:30)' * 0.2 /sqrt(20/(2 * 4));  %2.25 * (0:60)'/60;
% %rrr = 0.75 * 3 * (0:30)'/30;
% XXX = rrr*cos(ThetaXt);
% YYY = rrr*sin(ThetaXt);
% CCC = zeros(31,162);
% tempCCC = zeros(31,40);
% CCC(1:31,1:161) =  real( (spintimeBT20xi0p25ThetaB112Cut1SO3(1:31,:,4,2)+spintimeBT20xi0p25ThetaB112Cut1SO3(1:31,:,4,1)) ./ (spintimeBT20xi0p25ThetaB112Cut1SO3(1:31,:,1,2)+spintimeBT20xi0p25ThetaB112Cut1SO3(1:31,:,1,1)) );
% tempCCC = CCC(:,1:40);
% CCC(:,1:40) = CCC(:,81:120);
% CCC(:,81:120) = tempCCC;
% CCC(1:31,162) = CCC(1:31,2);
% contour(XXX,YYY,CCC, 20);
% hold on
% text(-4.75,4.4,str13, 'fontsize',14) 
% colorbar;
% axis equal
% axis([-5, 5,-5,5]);
% caxis([-0.5,0.5])
% xlhand = get(gca,'xlabel');
% ylhand = get(gca,'ylabel');
% set(xlhand,'string','$${x/\sqrt{Dt}}$$','Interpreter','latex','fontsize',14)
% set(ylhand,'string','$${y/\sqrt{Dt}}$$','Interpreter','latex','fontsize',14)
% set(gca,'XTick',[-5 0 5],'fontsize',14)
% set(gca,'YTick',[-5 0 5],'fontsize',14)
% set(gca,'XMinorTick','off')
% set(gca,'YMinorTick','off')
% set(gca,'ticklength',5*get(gca,'ticklength'))
% set(gca,'units','centimeters')
% set(gca,'Position',[2 2 3.8 3.8]);
% 
% 
% subplot(2,2,4);
% hold off 
% rrr = (0:60)' * 0.2 /sqrt(80/(2 * 4));  %2.25 * (0:60)'/60;
% XXX = rrr*cos(ThetaXt);
% YYY = rrr*sin(ThetaXt);
% CCC = zeros(61,162);
% tempCCC = zeros(61,40);
% CCC(1:61,1:161) =  real( (spintimeBT80xi0p25ThetaB112Cut1SO3(1:61,:,4,2)+spintimeBT80xi0p25ThetaB112Cut1SO3(1:61,:,4,1)) ./ (spintimeBT80xi0p25ThetaB112Cut1SO3(1:61,:,1,2)+spintimeBT80xi0p25ThetaB112Cut1SO3(1:61,:,1,1)) );
% tempCCC = CCC(:,1:40);
% CCC(:,1:40) = CCC(:,81:120);
% CCC(:,81:120) = tempCCC;
% CCC(1:61,162) = CCC(1:61,2);
% %CCC =  real( (spintimeBT80xi0p5ThetaB112Cut1SO2(1:81,:,4,2)+spintimeBT80xi0p5ThetaB112Cut1SO2(1:81,:,4,1)) ./ (spintimeBT80xi0p5ThetaB112Cut1SO2(1:81,:,1,2)+spintimeBT80xi0p5ThetaB112Cut1SO2(1:81,:,1,1)) );
% contour(XXX,YYY,CCC,20);
% hold on
% text(-4.75,4.4,str14, 'fontsize',14) 
% colorbar;
% axis equal
% axis([-5, 5,-5,5]);
% caxis([-0.05,0.05])
% xlhand = get(gca,'xlabel');
% ylhand = get(gca,'ylabel');
% set(xlhand,'string','$${x/\sqrt{Dt}}$$','Interpreter','latex','fontsize',14)
% set(gca,'XTick',[-5 0 5],'fontsize',14)
% set(gca,'YTick',[-5 0 5],'fontsize',14)
% set(gca,'YTickLabel',[])
% set(gca,'XMinorTick','off')
% set(gca,'YMinorTick','off')
% set(gca,'ticklength',5*get(gca,'ticklength'))
% set(gca,'units','centimeters')
% set(gca,'Position',[8.5 2 3.8 3.8]);
% 
% 
% set(gcf,'units','centimeters')
% set(gcf,'papersize',[5.8,4.6])
% set(gcf,'paperposition',[0,0,5.8,4.6])
% print('-depsc2', 'ProdPolarAllaCompareSO');
% 
% % str15 = sprintf('x10^{-3}');
% % text(8.9,5.2,str15, 'fontsize',14) 
% 
% NumThetasInPi2 = 20;
% LinearLength = 6;
% ThetaXt = zeros(1,(4 * NumThetasInPi2) + 2);
% for ttt=1:(4 * NumThetasInPi2) + 2 
%        ThetaXt(ttt) = ((ttt-1) * 2 * pi)/(4 * NumThetasInPi2);
% end
% XXX = zeros(LinearLength,(4 * NumThetasInPi2) + 2);
% YYY = zeros(LinearLength,(4 * NumThetasInPi2) + 2);
% CCC = zeros(LinearLength,(4 * NumThetasInPi2) + 2);
% rrr = (0:LinearLength-1)';  %2.25 * (0:60)'/60;
% XXX = rrr*cos(ThetaXt);
% YYY = rrr*sin(ThetaXt);
% % % Sx to Sz, I think.
% % CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polarE40B0Tau15(4,2,1:LinearLength,:)+polarE40B0Tau15(4,1,1:LinearLength,:)) ./ (polarE40B0Tau15(1,2,1:LinearLength,:)+polarE40B0Tau15(1,1,1:LinearLength,:)) );
% % Sx to Sy, I think.
% %  CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polar6repeat(3,2,1:LinearLength,:)+polar6repeat(3,1,1:LinearLength,:)) ./ (polar6repeat(1,2,1:LinearLength,:)+polar6repeat(1,1,1:LinearLength,:)) );
% % % %I to Sy, I think.
% % CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polarE40B0Tau15(3,1,1:LinearLength,:)) ./ (polarE40B0Tau15(1,1,1:LinearLength,:)) );
% % %I to Sz, I think.
% % CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =  real( (polarE40B0Tau15(4,1,1:LinearLength,:)) ./ (polarE40B0Tau15(1,1,1:LinearLength,:)) );
% % % %I to I, I think.
% CCC(1:LinearLength,1:(4 * NumThetasInPi2) + 1) =     (polar6red(1,1,1:LinearLength,:));
% tempCCC = CCC(:,1:NumThetasInPi2);
% CCC(:,1:NumThetasInPi2) = CCC(:,(2*NumThetasInPi2)+1:3*NumThetasInPi2);
% CCC(:,(2*NumThetasInPi2)+1:3*NumThetasInPi2) = tempCCC;
% CCC(:,(4*NumThetasInPi2)+2) = CCC(:,2);
% CCC(:,(4 * NumThetasInPi2) + 2)=CCC(:,2);
% CCC(:,(4 * NumThetasInPi2) + 1)=CCC(:,1);
% contour(XXX,YYY,max(abs(CCC),12),20);
% colorbar;
