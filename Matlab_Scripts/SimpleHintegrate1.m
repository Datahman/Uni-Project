function [DOS, kfAtThetaZero, KfAtThetaPi2, RhoN, DinverseTot] =  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumGridPoint,MaxNumEf,NumValueTheta, tau)

%variables calculated by the loop over thetak
thetaarray = zeros(NumValueTheta -1,1 ); % values of the angle theta
DOSAtThetaArray = zeros(NumValueTheta -1,1 ); %at each angle thetak, Scalar version of the DOS
%DOSAtThetaArrayA = zeros(NumValueTheta -1,1);
%DOSAtThetabArray = zeros(NumValueTheta -1, 1);
%KfIndArray = zeros(NumValueTheta -1, 2); %  for two eigen-values
%DerivtiveArray = zeros(NumValueTheta -1, 2); % Array for derivatives
RhoN = zeros(Norbitals,Norbitals); % Matrix density of states.
RhoTmp = zeros(Norbitals,Norbitals); % Matrix density of states.
% Diffusion operator, at each value of thetak
D_inverse = zeros(Norbitals * Norbitals, Norbitals * Norbitals,NumValueTheta -1); % Diffuson operator 16x16
D_inverse1 = zeros(Norbitals * Norbitals, Norbitals * Norbitals); % Diffuson operator 16x16
        n1 = zeros(NumValueTheta -1,1 );   % used for plotting 
n2 = zeros(NumValueTheta -1,1 );
n3 = zeros(NumValueTheta -1,1 );
n4 = zeros(NumValueTheta -1,1 );
n5 = zeros(NumValueTheta -1,1 );
kfs = zeros(NumValueTheta - 1,MaxNumEf); % stores kfs at each value of theta


% variables used inside the loop over thetak
FSNumArray = zeros(NumValueTheta -1,1 ); % at each angle thetak, the number of Fermi surfaces
kfArray1 = zeros(MaxNumEf, NumValueTheta-1); % at each value of thetak, the kfs and ?

kfArray = zeros(MaxNumEf,3); % at each value of thetak, the kfs and ?

% ProbVector is the vector that for the total probability
ProbVector = zeros(Norbitals * Norbitals, 1);
for nn = 1:Norbitals
    ProbVector(nn + (Norbitals * (nn - 1))) = 1;
end

for n = 1: NumValueTheta - 1  % loop over the angle around the origin
    thetak = (2* pi / (NumValueTheta - 1)) * (n-1); % decide on the angle thetak
    thetaarray(n) = thetak; % store the angle in an array

    % get all the Fermi surfaces, and a count of them
    % todo Bug Bug Bug: we are getting wierd fermi surfaces including
    % negative values.
    [FermiSurfaceNum,kfArray] = excitonFermiSurfaceA(Norbitals, thetak,BMagnetic,deltab,Ef,kmax,NumGridPoint,MaxNumEf) ;  

 
        countcont = 0; % check that we are getting the right Fermi surfaces

    % These are for analytical calculation of the Fermi surface in the case
    % of a quadratic energy
%     DOSAtThetaAnalytic =0;
%     kfd1 = - BMagnetic * sin(thetak) + sqrt((BMagnetic * BMagnetic) * (sin(thetak) * sin(thetak)) + (Ef - BMagnetic* BMagnetic));
%     kfd2 = - BMagnetic * sin(thetak) - sqrt((BMagnetic * BMagnetic) * (sin(thetak) * sin(thetak)) + (Ef - BMagnetic* BMagnetic));
%     KfIndArray(n,1) = kfd1;
%     KfIndArray(n,2) = kfd2;
%     AnalyticalDerivative1 = 2 * kfd1 + 2 * BMagnetic * sin(thetak); % dE_{i} / dk_f , i = 1,2 
%     AnalyticalDerivative2 = 2 * kfd2 + 2 * BMagnetic * sin(thetak);
%     DerivtiveArray(n,1) = AnalyticalDerivative1;
%     DerivtiveArray(n,2) = AnalyticalDerivative2;
    
    DOSAtTheta =0; % Initialize the scalar version of the DOS
%    DOSAtThetab =0; % check variable
    
    for iFS=1:FermiSurfaceNum % sum over the Fermi surfaces
 
        % find this Fermi surface's contribution to the scalar density of
        % states
    DOSAtTheta = DOSAtTheta + (( (2 * pi)/ (NumValueTheta -1)) * (kfArray( iFS,1) / abs((kfArray(iFS,3)))));

    % get the Hamiltonian and derivatives at this Fermi surface
    hhha = excitonHamiltonian(kfArray(iFS,1), thetak, BMagnetic, deltab);
    % get the eigenvalues and eigenvectors of the Hamiltonian at this Fermi surface
    [vv,dd] = eig(hhha(:,:,1));

    
    for r=1:Norbitals % look for the eigenvalue that matches the Fermi surface that we are on
        if abs(Ef - dd(r,r)) < 0.00000001  % if statement for matching the Fermi surface
            countcont = countcont + 1;
            vvt = vv';% transpose
           
            % calculate the matrix version of the density matrix
            RhoTmp = (( (2 * pi)/(NumValueTheta-1) )*((vv(:,r) * vvt(r,:) ) * (kfArray( iFS,1) / abs((kfArray(iFS,3))))));
            RhoN = RhoN + RhoTmp;
            %trace(vv(:,m) * vvt(m,:))
            %DOSAtThetabArray(n) = DOSAtThetabArray(n)+(trace((vv(:,r) * vvt(r,:) ) * (kfArray( iFS,1) / abs((kfArray(iFS,3))))));

            % todotodo: 
            % 2. Restructure loops so that the r loop is inside the
            % s1,s2,s3,s4 loops.
            % 3. Add q-dependence (momentum dependence) and produce plots
            % of the exciton distribution similar to the plots in the
            % polariton paper.
 D_inverse1 = zeros(Norbitals * Norbitals, Norbitals * Norbitals);
 
            % Four loops over all  matrix elements in the diffuson
            for s1 = 1:Norbitals
    for s2 = 1:Norbitals
        for s3 = 1:Norbitals
            for s4 = 1:Norbitals

                t1 = (s1 - 1) * 4 + s2;
                t2 = (s3 - 1) * 4 + s4;
                    for m = 1:Norbitals
                        delta_plus = dd(m,m) - dd(r,r);
                        D_inverse1(t1,t2) = D_inverse1(t1,t2) + (( (1 * pi)/(NumValueTheta-1) ) * vvt(r,s1) * vv(s3,r) * vv(s2,m) * vvt(m,s4) * (kfArray( iFS,1) / abs((kfArray(iFS,3))))) / (1- 1i * tau * delta_plus);
                        D_inverse1(t1,t2) = D_inverse1(t1,t2) + (( (1 * pi)/(NumValueTheta-1) ) * vvt(m,s1) * vv(s3,m) * vv(s2,r) * vvt(r,s4) * (kfArray( iFS,1) / abs((kfArray(iFS,3))))) / (1 + 1i * tau * delta_plus);
                     end

                    % end of four loops over all matrix elements in the diffuson
            end
        end
    end
            end
             D_inverse(:,:,n) = D_inverse(:,:,n) + D_inverse1;

%        check whether multiplying by ProbVector gives Rho
 res = zeros(Norbitals * Norbitals);
 resm = zeros(Norbitals, Norbitals);
            res = D_inverse1 * ProbVector;
            for s1 = 1:Norbitals
                for s2 = 1:Norbitals
             t1 = (s2 - 1) * 4 + s1;
                    resm(s1, s2) = res(t1);
                end
            end
            if norm(resm-RhoTmp) > 0.00000001
                 norm(res), resm, RhoTmp, resm-RhoTmp
                 %                D_inverse1
                %ProbVector
                msgbox "probability not conserved."
            end
            
            
        end  % end of if statement for matching the Fermi surface
    end  % end of loop that looks for the eigenvalue that matches the Fermi surface
    
        
    end % end of loop over Fermi surfaces

    % check how many Fermi surfaces were found
    if countcont < FermiSurfaceNum
        countcont
        FermiSurfaceNum
%         msgbox "wrong count of fermi surfaces"
    end

    % if(abs(DOSAtThetabArray(n) - DOSAtTheta) > 0.0001)
   % t
   % n
  %  kfArray
 %   DOSAtThetabArray(n)
%    DOSAtTheta
    
    %DOSAtTheta
    %end
%    for s = 1:2
       % DOSAtThetaAnalytic = DOSAtThetaAnalytic +  abs(((KfIndArray(n,s))/DerivtiveArray(n,s)));
%    end

        %DOSAtThetaArrayA(n) =  DOSAtThetaAnalytic;
        DOSAtThetaArray(n)  = DOSAtTheta;   

        % these statements prep for plotting
kfArray1 (:,n) = kfArray(:,1); 

n1(n)=thetaarray(n);
n2(n) = kfArray1(1,n); % First eig value 
n3(n) = kfArray1(2,n);% second eig value
n4(n) = kfArray1(3,n);% second eig value
n5(n) = kfArray1(4,n);% second eig value
%kfs(n) = kfArray1(3,3,n); % third eigstate
%n5(n) = kfArray1(4,3,n); %4th eigstate

% size(kfArray(:))
% size(kfArray1(:,3,n))


end % loop over theta
%size(kfArray )
%size(DerivtiveArray) 


% this is good graphing code for the Fermi surfaces
% hold off;
% polar(n1,n2 , 'red');
% hold on;
% polar(n1,n3, 'yellow');
% polar(n1,n4, 'yellow');
% hold on;
% polar(n1,n5, 'black');



%hold on;
%kfArray;
%DOSAtThetaArrayA  = DOSAtThetaArrayA  ;
%polar(n1,kfs(:,1), 'blue');
%hold on;
%polar(n1,kfs(:,2), 'red');
%hold on;
%polar(n1,kfs(:,3), 'black');
%hold on;
%polar(n1,kfs(:,4), 'yellow');

% polar(n1,DOSAtThetaArray, 'red');

TotalDOSatE=sum(DOSAtThetaArray) ;
%TotalDOSatEAnalytic = (sum(DOSAtThetaArrayA) * (2 * pi)/ (NumValueTheta -1)); 
%TotalDOSatE2=sum(DOSAtThetabArray) * (2 * pi)/ (NumValueTheta -1);
%DOSAtThetaArray - DOSAtThetabArray
%TotalDOSatE2 - TotalDOSatE
DOS = TotalDOSatE;

% DOS
% 
% trace(RhoN)
%DOSA = TotalDOSatEAnalytic; % Returns 1x1 .. to be fixed.
%n1lim = 2;
%axis([-1 1 -1 1] * n1lim);

DinverseTot = zeros(Norbitals * Norbitals, Norbitals * Norbitals); % Diffuson operator 16x16
for n = 1: NumValueTheta - 1  % loop over the angle around the origin
 DinverseTot = DinverseTot + D_inverse(:,:,n) ;
end
% res1 = zeros(Norbitals * Norbitals);
% res1res = zeros(Norbitals * Norbitals);
% res1 = DinverseTot * ProbVector;
% res1res = res1 - (ProbVector' * res1) * ProbVector / (ProbVector' * ProbVector);
% if norm(res1res)
%     DinverseTot
%     res1, res1res
%     msgbox 'Probability not conserved 2.'
% end

%        check whether multiplying by ProbVector gives Rho
 res = zeros(Norbitals * Norbitals);
 resm = zeros(Norbitals, Norbitals);
            res = DinverseTot * ProbVector;
            for s1 = 1:Norbitals
                for s2 = 1:Norbitals
             t1 = (s2 - 1) * 4 + s1;
                    resm(s1, s2) = res(t1);
                end
            end
            if norm(resm-RhoN) > 0.00000001
                 norm(res), resm, RhoN, resm-RhoN
                 %                D_inverse1
                %ProbVector
                msgbox "probability not conserved 2."
            end
%             trace(RhoN)
% trace(RhoN)/Norbitals
% todo consider whether to divide by RhoN rather than by the trace of RhoN.
%  Maybe do the division in a symmetrized way.
            DinverseTot = DinverseTot / (trace(RhoN)/Norbitals);
            DinverseTot = eye(Norbitals * Norbitals) - DinverseTot;
            
kfAtThetaZero = kfArray1(:, 1); % For thetea = 0
KfAtThetaPi2  = kfArray1(:,floor((NumValueTheta/ 4))); % For theta = approx pi/2.
%KfAtThetaPi2
end