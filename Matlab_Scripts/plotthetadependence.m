

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

%plotthetadependence(4, 0, 0.5, 40, 100, 100, 100, 8, 1)

function [abc] = plotthetadependence(Norbitals, BMagnetic,deltab, Ef, NumThetas, kmax,NumKGridPoints, MaxNumFS, PlotType )

%todo: also get the rho eigenvalues

ThetaArray = zeros(NumThetas,1); 

% Ef should be calculated relative to Ef = -10 because the dispersion's
% minimum is around -10.
FSNumArray = zeros(NumThetas,1);
kfArray = zeros(NumThetas,MaxNumFS); % 1/ micro m
% vf / hbarOverTau determines the scattering length
vfArray = zeros(NumThetas,MaxNumFS); % micro eV  * micro m
% deltaE / hbarOverTau determines EY vs. DP
% deltaE / Ef determines whether SO is strong
% Ef / hbarOverTau must always be much larger than 1 to prevent localization physics.   
deltaEArray = zeros(NumThetas, MaxNumFS,Norbitals- 1); % micro eV
% The spin precession length is vf / deltaE
lSOArray = zeros(NumThetas, MaxNumFS, Norbitals-1);

DegeneracyArray = zeros(NumThetas,MaxNumFS);


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

for ThetaCounter= 1:NumThetas
    ThetaArray(ThetaCounter) =  ((2 * pi)* (ThetaCounter -1 )/(NumThetas - 1)); 
 [FSNum,kfArray1] = excitonFermiSurfaceA(Norbitals, ThetaArray(ThetaCounter),BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumFS) ;  
   FSNumArray(ThetaCounter) = FSNum;
   kfArray(ThetaCounter,:) = kfArray1(:,1);
   vfArray(ThetaCounter,:) = -1;
   DegeneracyArray(ThetaCounter,:) = kfArray1(:,6);
   
   for FSCounter=1:FSNum

       if DegeneracyArray(ThetaCounter,FSCounter) < 0
           continue;
       end

       
       % get the Hamiltonian and velocity operators
  hhh1 = excitonHamiltonian(kfArray(ThetaCounter,FSCounter),ThetaArray(ThetaCounter), BMagnetic, deltab);
  % get the eigenvalues
[vv,dd] = eig(hhh1(:,:,1));

            % This calculates the velocity operator along the angle ThetaArray(ThetaCounter).
            hhh2 = (hhh1(:,:,2) * cos(ThetaArray(ThetaCounter))) + (hhh1(:,:,3) * sin(ThetaArray(ThetaCounter)));
        vvt = vv';         
        vfArray1 = zeros(Norbitals);
        % compute the Fermi velocities
        for EigCounter = 1:Norbitals
         vfArray1(EigCounter) = vvt(EigCounter,:) * hhh2 *  vv(:, EigCounter);
        end
        
        % find the eigenvalue for this FS
        iWhichFS = 0;
        for EigCounter = 1:Norbitals
            if( abs(dd(EigCounter,EigCounter) - Ef) < 1e-13)
                iWhichFS = EigCounter;
                break %leave the loop since we found the FS
            end
        end
        if iWhichFS < 1
            msgbox 'Did not find any eigenvalue for this FS'
        end
            
         % todo this does not take into account that there may be several
         % degenerate states, each with a different vf.
         vfArray(ThetaCounter,FSCounter) = vfArray1(iWhichFS);

         % compute the deltaE's
        EigCounterAdjustment = 0;
        for EigCounter = 1:Norbitals
            if EigCounter == iWhichFS
                EigCounterAdjustment = -1;
                continue % skip this case
            end
            deltaEArray(ThetaCounter, FSCounter,EigCounter+EigCounterAdjustment) = dd(EigCounter,EigCounter) - Ef; 
            lSOArray(ThetaCounter, FSCounter,EigCounter+EigCounterAdjustment) = 1/((dd(EigCounter,EigCounter) - Ef)/vfArray1(EigCounter)); 
        end
   end
end

switch PlotType
    case 1
% this plots the kfs
hold off;
for FSCounter = 1:MaxNumFS
    polar(ThetaArray(:,1),kfArray(:,FSCounter));
    xlabel('\theta_{k}')
    ylabel('$\vec{k}_{F}$','interpreter', 'latex')
    hold on;
end

    case 2
% this plots the Fermi velocities
hold off
for FSCounter = 1:MaxNumFS
    polar(ThetaArray(:,1),real(vfArray(:,FSCounter)),'r');
    xlabel('\theta_{k}')
    ylabel('$ \vec{v}_{F}$','interpreter', 'latex')
    hold on;
end
   
    case 3
% this plots the energy splittings
hold off;
for FSCounter = 1:MaxNumFS
    for OrbCounter = 1:Norbitals-1
        plot(ThetaArray(:,1),deltaEArray(:,FSCounter,OrbCounter));
         xlabel('\theta_{k}')
         ylabel('$ \Delta E$','interpreter', 'latex')
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
        plot(ThetaArray(:,1),abs(lSOArray(:,FSCounter,OrbCounter)), 'r');
        hold on;
    end
end
    axis([0, 2 * pi,0,2]);

    case 5
% this plots the Fermi surface numbers
hold off;
    polar(ThetaArray(:),FSNumArray(:));
    xlabel('\theta_{k}')
    ylabel('FSN')
    hold on;
    otherwise
        msgbox 'Bad plottype';
end

end