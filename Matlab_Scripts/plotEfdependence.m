% % first run with 20 pts in the angular integration - previously had 8 and generated clear octopole patterns.  Still scars around r = 9, but they moved close to the y axis and the whole graph looks a bit more controlled.  Still a 20-pole pattern, and very fast oscillations along the radial axis.
% polarE40B0Tau15HiResB = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 40, 20,4);
% % 20 pts in the angular integration, but a rougher lattice.  This gives results that are wildly different, and less singular.
% polarE40B0Tau15HiResC = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.3, 20, 20,4);
% % This should be a quarter of HiResB but there is no clear connection.
% polarE40B0Tau15HiResD = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 20, 20,4);
% % shows the original square grid, which oscillates at every other lattice unit
% [polarE40B0Tau15HiResD,spinxyHiResD] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 20, 1,4);
% % does not do the final FFT - found that the 1,1 point is replicated to 11,1 and 1,11 and 11,11
% [polarE40B0Tau15HiResE,spinxyHiResE] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 20, 1,4);
% % does not do the exponential or the FFT - the problem is that both the cos and the sin are set to zero at these points, even though the magnitude is not zero.
% [polarE40B0Tau15HiResF,spinxyHiResF] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 20, 1,4);
% % now trying with odd numbers, maybe the problem will go away, intensity tightly concentrated in the center, 4800 seconds
% [polarE40B0Tau15HiResG,spinxyHiResG] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,4);
% % seeing whether it spreads out when time is increased -yes, here it is flat, good
% [polarE40B0Tau15HiResH,spinxyHiResH] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,36);
% % checking convergence with lattice spacing - almost identical to H, good
% [polarE40B0Tau15HiResI,spinxyHiResI] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.075, 39, 20,36);
% % making a time series of diffusion G-K-J-L-H. t=16,24,36 are almost flat, i.e. completely dominated by finite size effects caused by the small sample
% [polarE40B0Tau15HiResK,spinxyHiResK] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,8);
% [polarE40B0Tau15HiResJ,spinxyHiResJ] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,16);
% [polarE40B0Tau15HiResL,spinxyHiResL] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,24);
% % checking convergence with lattice spacing again - Sx to Sz and I to Sz is stable, the I to I I'm not sure about, and the Sx to Sy and I to Sy is not stable
% [polarE40B0Tau15HiResM,spinxyHiResM] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.075, 39, 20,4);
% [polarE40B0Tau15HiResN,spinxyHiResN] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.075, 39, 20,8);
% % extending the time series to shorter times O-P-G-K
% [polarE40B0Tau15HiResO,spinxyHiResO] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,1);
% [polarE40B0Tau15HiResP,spinxyHiResP] = IXPolarPlot(4,0,0.5,40, 100,100,8,20,15, 0.15, 19, 20,2);
% % extending the time series to B=7 Q-R-S-T
% [polarE40B0Tau15HiResQ,spinxyHiResQ] = IXPolarPlot(4,7,0.5,40, 100,100,8,20,15, 0.15, 19, 20,1);
% [polarE40B0Tau15HiResR,spinxyHiResR] = IXPolarPlot(4,7,0.5,40, 100,100,8,20,15, 0.15, 19, 20,2);
% [polarE40B0Tau15HiResS,spinxyHiResS] = IXPolarPlot(4,7,0.5,40, 100,100,8,20,15, 0.15, 19, 20,4);
% [polarE40B0Tau15HiResT,spinxyHiResT] = IXPolarPlot(4,7,0.5,40, 100,100,8,20,15, 0.15, 19, 20,8);
% [polarE40B0Tau15HiResU,spinxyHiResU] = IXPolarPlot(4,7,0.5,40, 100,100,8,20,15, 0.075, 39, 20,4);


% todo: what is the Ioffe-Regel criterion, exactly? 

% polarE40B0Tau15a = IXPolarPlot(4,0,0.5,40, 100,100,8,8,15, 0.6, 18, 8,4);

% polarE40B0Tau15 = IXPolarPlot(4,0,0.5,40, 100,100,8,8,15, 0.6, 12, 8,4);

% polarE40B0Tau15 = IXPolarPlot(4,0,0.5,40, 100,100,8,8,15, 0.3, 40, 20,10);
% polarE1000B0Tau300 = IXPolarPlot(4,0,0.5,1000, 100,100,8,8,300, 0.07, 40,20, 10);

%polarE40B0Tau15 = IXPolarPlot(4,0,0.5,40, 100,100,8,8,15, 0.3, 40, 20,4);
%polarE40B0Tau15fast = IXPolarPlot(4,0,0.5,40, 100,100,8,8,15, 1.5, 8, 20,4);
% At E = 40, there are deltaE's ranging from 30, 50, 80, ...200
% Choose hbarOverTau = 15
% At E = 40, all the vfs are about the same: vf=6
% At E = 40, lSO ranges from 0.02 up to 0.25, and filling up that range
% vf / hbarOverTau = scattering length = 0.4 
% The linear length of the sample must be well over 4l = 1.6, multiplied by
% the square root of t/tau
% deltaE / hbarOverTau = 2 ... 7 , so pretty much in the EY regime
% The spin decay time is the same as Tau, within factors of O(1).  
% The spin precession length is lSO = vf / deltaE 
%                                = (vf / hbarOverTau) / (deltaE / hbarOverTau)
% The spin decay length is the same as the scattering length, 0.4 
% the resolution should be 0.15
% the linearlength should be well over (1.6/0.15 ) * sqrt(t/tau) 
% t should be a few multiples of Tau

% At E = 500, there are deltaE's at 100, 200, 300, 400.  
% Choose hbarOverTau = 150
% At E = 500, all the vfs are about the same: vf=20
% At E = 500, lSO ranges from 0.03 up to 0.22
% vf / hbarOverTau = scattering length = 2/15 = 0.13 
% deltaE / hbarOverTau = 0.7 ... 2.7, so near the EY-DP transition, on the
% EY side.
% The spin decay time is the same as Tau, within factors of O(1-5).  
% The spin decay length is the same as or a bit bigger than the scattering
% length, 0.13 - 0.5
% the resolution should be 0.05, and the radius should be 1 or 5 or so, and
% t should be a few to tens of  Tau


%polarE1000B0Tau300 = IXPolarPlot(4,0,0.5,1000, 100,100,8,8,300, 0.07, 40,20, 4);
%polarE1000B0Tau300 = IXPolarPlot(4,0,0.5,1000, 100,100,8,8,300, 0.35, 8, 20, 4);
% At E = 1000, there are deltaE's at 150, 300, 450, 600.  
% Choose hbarOverTau = 300
% At E = 1000, all the vfs are about the same: vf=27
% At E = 1000, lSO ranges from 0.04-0.12, and then 0.2
% vf / hbarOverTau = scattering length = 27/300 = 0.09 
% The linear length of the sample must be well over 4l = 0.36, multiplied by
% the square root of t/tau
% deltaE / hbarOverTau = 0.5 ... 2, so near the EY-DP transition
% The spin decay time is the same as Tau, within factors of O(1-5).  
% The spin decay length is the same as or a bit bigger than the scattering
% length, 0.09 - 0.5
% the resolution should be 0.05
% the linearlength should be well over (0.36/0.05 ) * sqrt(t/tau) 
% t should be a few to tens of  Tau

%plotEfdependence(4, 0,0, 0.5, 10, 100, 100, 100, 100, 8,1)

function [abc] = plotEfdependence(Norbitals, thetak,BMagnetic,deltab, Efmin, Efmax, NumEfs, kmax,NumKGridPoints, MaxNumFS, PlotType )

%todo: also get the rho eigenvalues

EfArray = zeros(NumEfs,1); 

% Ef should be calculated relative to Ef = -10 because the dispersion's
% minimum is around -10.
FSNumArray = zeros(NumEfs,1);
kfArray = zeros(NumEfs,MaxNumFS); % 1/ micro m
% vf / hbarOverTau determines the scattering length
vfArray = zeros(NumEfs,MaxNumFS); % micro eV  * micro m
% deltaE / hbarOverTau determines EY vs. DP
% deltaE / Ef determines whether SO is strong
% Ef / hbarOverTau must always be much larger than 1 to prevent localization physics.   
deltaEArray = zeros(NumEfs, MaxNumFS,Norbitals- 1); % micro eV
% The spin precession length is vf / deltaE
lSOArray = zeros(NumEfs, MaxNumFS, Norbitals-1);

DegeneracyArray = zeros(NumEfs,MaxNumFS);

% In the EY regime deltaE / hbarOverTau >> 1:
% The spin decay time is the same as Tau, within factors of O(1).  
% The spin precession length is vf / deltaE 
%                                = (vf / hbarOverTau) / (deltaE / hbarOverTau)
% The spin decay length is the same as the scattering length, vf / hbarOverTau. 

% In the DP regime deltaE / hbarOverTau << 1:
% The spin decay time is the same as Tau / (deltaE / hbarOverTau)^2 , within factors of O(1).  
% The spin precession length is vf / deltaE 
%                                = (vf / hbarOverTau) / (deltaE / hbarOverTau)
% The spin decay length is the same as the spin precession length. 

kfArray1 = zeros(MaxNumFS,6);
 hhh1 = zeros(Norbitals, Norbitals,3);

for ECounter= 1:NumEfs
    EfArray(ECounter) = Efmin + ((Efmax - Efmin)* (ECounter -1 )/(NumEfs - 1)); 
 [FSNum,kfArray1] = excitonFermiSurfaceA(Norbitals, thetak,BMagnetic,deltab,EfArray(ECounter),kmax,NumKGridPoints,MaxNumFS) ;  
   FSNumArray(ECounter) = FSNum;
   kfArray(ECounter,:) = kfArray1(:,1);
   vfArray(ECounter,:) = -1;
   DegeneracyArray(ECounter,:) = kfArray1(:,6);
%    FSNum
%    kfArray(ECounter,:)
   
   for FSCounter=1:FSNum
       
       if DegeneracyArray(ECounter,FSCounter) < 0
           continue;
       end
           % get the Hamiltonian and velocity operators
  hhh1 = excitonHamiltonian(kfArray(ECounter,FSCounter), thetak, BMagnetic, deltab);
  %hhh1 = simpleHamiltonian(kfArray(ECounter,FSCounter), thetak, BMagnetic, deltab);

  % get the eigenvalues
[vv,dd] = eig(hhh1(:,:,1));

            % This calculates the velocity operator along the angle thetak.
            hhh2 = (hhh1(:,:,2) * cos(thetak)) + (hhh1(:,:,3) * sin(thetak));
        vvt = vv';         
        
        vfArray1 = zeros(Norbitals);
        % compute the Fermi velocities
        for EigCounter = 1:Norbitals
         vfArray1(EigCounter) = vvt(EigCounter,:) * hhh2 *  vv(:, EigCounter);
        end
        
        % find the eigenvalue for this FS
        iWhichFS = 0;
        for EigCounter = 1:Norbitals
            if( abs(dd(EigCounter,EigCounter) - EfArray(ECounter)) < 1e-13)
                iWhichFS = EigCounter;
                break %leave the loop since we found the FS
            end
        end
         if iWhichFS < 1
            msgbox 'Did not find any eigenvalue for this FS'
        end
       
         % todo this does not take into account that there may be several
         % degenerate states, each with a different vf.
         vfArray(ECounter,FSCounter) = vfArray1(iWhichFS);

        % compute the deltaE's
        EigCounterAdjustment = 0;
        for EigCounter = 1:Norbitals
            if EigCounter == iWhichFS
                EigCounterAdjustment = -1;
                continue % skip this case
            end
            deltaEArray(ECounter, FSCounter,EigCounter+EigCounterAdjustment) = dd(EigCounter,EigCounter) - EfArray(ECounter); 
            lSOArray(ECounter, FSCounter,EigCounter+EigCounterAdjustment) = 1/((dd(EigCounter,EigCounter) - EfArray(ECounter))/vfArray1(EigCounter)); 
        end
   end
end

switch PlotType
    case 1
% this plots the kfs
hold off;
for FSCounter = 1:MaxNumFS
    
    plot(EfArray(:,1),kfArray(:,FSCounter));
     xlabel('Energy ($\mu eV$)','interpreter', 'latex')
     ylabel('$\vec{k_{F}}$','interpreter', 'latex')
    hold on;
   
end
    case 2
%this plots the Fermi velocities
hold off
for FSCounter = 1:MaxNumFS
    plot(EfArray(:,1),vfArray(:,FSCounter),'r');
     xlabel('Energy ($\mu eV$)','interpreter', 'latex')
     ylabel('$\vec{v_{F}}$','interpreter', 'latex')
    hold on;
end

    case 3
%this plots the energy splittings
hold off;
for FSCounter = 1:MaxNumFS
    for OrbCounter = 1:Norbitals-1
        plot(EfArray(:,1),deltaEArray(:,FSCounter,OrbCounter));
          xlabel('Energy ($\mu$ eV)','interpreter', 'latex')
          ylabel('$\Delta$ E','interpreter', 'latex')
        hold on;
    end
end

    case 4
%this plots the spin-orbit decay lengths, I think.
% todo when we calculate this we get some infinites because of dividing by
% zero.  We should fix this before publishing any graph based on this.
hold off;
for FSCounter = 1:MaxNumFS
    for OrbCounter = 1:Norbitals-1
        plot(EfArray(:,1),abs(lSOArray(:,FSCounter,OrbCounter)), 'r');
        hold on;
    end
    axis([-20,100,0,2]);
end
    case 5
% this plots the Fermi surface numbers
hold off;
    plot(EfArray(:),FSNumArray(:));
      xlabel('Energy ($\mu$ eV)','interpreter', 'latex')
      ylabel('FSN','interpreter', 'latex')
    hold on;

    otherwise
        msgbox 'Bad plottype';
end
end