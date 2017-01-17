
% [IntegratedDOSa, kfAtThetaZeroa, KfAtThetaPi2a, IntegratedRhoNa, IntegratedDinversea, DiffExpCoeffsa] =  SimpleHintegrate(4,0,0.5,40,100,100,8,10, 15, kspacegrid, 1);


% this function calculates the scalar and matrix DOS, and Dinverse, after integration
% over the Fermi surface.
% It also calculates the Fermi momenta at theta = 0 and theta = pi/2. 
% It does this at particular values of Ef, tau, and the Hamiltonian parameters.  
% todo: someday add code that handles the sqrt divergence that occurs when
% a Fermi surface does not encircle the origin.
function [IntegratedDOS, kfAtThetaZero, KfAtThetaPi2, IntegratedRhoN, IntegratedDinverse, DiffExpCoeffs] =  SimpleHintegrate(Norbitals,BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumEf,NumValueTheta, hbarOverTau, kspacegrid, LinearLength)
 kspacegrid=zeros(1,1,3);
% % the following block is for checking the validity of our approximations.
% me =  0.511 * 10^12; % mass of the electron, micro eV
% mex = 0.21 * me;  % exciton mass, in rough agreement with Szymanska and Littlewood
% hbarc = 1.24 * 10^6 / (2 * pi); %micro eV * micro m
% hbarvf = hbarc * (hbarc * kf) / mex;  % the Fermi velocity if the spin-orbit coupling is neglected
% Ef = (hbarc * kf)^2 / (2 * mex); % the Energy if the spin-orbit coupling is neglected. 43.0 for Butov's numbers.  
% % the binding energy of the exciton is in the range of 4000 micro eV from
% % Szymanska and Littlewood, so Ef should be less than 1300 micro eV or so.
% Eprcsn = splittingbetweenenergylevels; %micro eV % should be the splitting between the max and min eigenvalues
% EprcsnEf = Eprcsn/Ef; % should be small, otherwise can expect spin polarization, ??? perhaps only if B != 0 ???
% % With Butov's numbers the spin-orbit eigenvalues (without p^2/2m, with B=7T) are
% % +-55.9, and +- 27.4, and don't depend on thetak.  These numbers changes by < 0.2 if B=0T.
% lll = hbarvf  / ETau; % this is the scattering length, tau = hbar / ETau
% xi = Eprcsn / ETau; % small in the DP regime, large in the EY regime.
% alpha = ETau / Ef; % should be small, otherwise becomes localized and the diffuson is wrong

IntegratedDOS = 0;
IntegratedRhoN = zeros(Norbitals, Norbitals);
IntegratedCollisionOperator = zeros(LinearLength,LinearLength,Norbitals * Norbitals, Norbitals * Norbitals);
NumDiffExpCoeffs = 6;
DiffExpCoeffs = zeros(NumDiffExpCoeffs, Norbitals * Norbitals, Norbitals * Norbitals);

%array variables calculated by the loop over thetak
thetaArray = zeros(1,NumValueTheta); % values of the angle theta
 kfsArray = zeros(MaxNumEf,NumValueTheta); % stores kfs at each value of theta
FSNumArray = zeros(1,NumValueTheta); % at each angle thetak, the number of Fermi surfaces
DOSAtThetaArray = zeros(1,NumValueTheta); %scalar version of the DOS
RhoNArray = zeros(Norbitals, Norbitals, NumValueTheta);

RhoTmp = zeros(Norbitals,Norbitals); % Matrix density of states - working variable.
 % the Fermi momentum and x and y components of the Fermi velocity
kfArray = zeros(MaxNumEf,6);




% ProbVector is the vector for the total probability
ProbVector = zeros(Norbitals * Norbitals, 1);
for nn = 1:Norbitals
    ProbVector(nn + (Norbitals * (nn - 1))) = 1;
end

for n = 1: NumValueTheta  % loop over the angle around the origin
    thetak = (2* pi / (NumValueTheta)) * (n-1); % decide on the angle thetak
    thetaArray(n) = thetak; % store the angle in an array

%     thetak
    % get all the Fermi surfaces, the associated Fermi velocities, and a
    % count of the Fermi surfaces
% This function returns the count of Fermi surfaces in FermiSurfaceNum.
% kfArray's first element is kf, its second element is an index into the
% array of k's, and its third element is the Fermi velocity at this kf and
% thetak.
    [FSNumArray(n),kfArray] = excitonFermiSurfaceA(Norbitals, thetak,BMagnetic,deltab,Ef,kmax,NumKGridPoints,MaxNumEf) ;  
kfsArray(:,n) = kfArray(:,1);

    
    for iFS=1:FSNumArray(n) % sum over the Fermi surfaces
        if kfArray(iFS,6) < 0 % this is for not handling a degeneracy twice
            continue;
        end
        
        
    % get the Hamiltonian and derivatives at this Fermi surface
    hhha = excitonHamiltonian(kfArray(iFS,1), thetak, BMagnetic, deltab);
    % get the eigenvalues and eigenvectors of the Hamiltonian at this Fermi surface
    [vv,dd] = eig(hhha(:,:,1));
    
        NumDegenerateStates=0;  
        DegenerateStateIndices = zeros(MaxNumEf);
        for ll=1:Norbitals
            if abs(Ef - dd(ll,ll)) < 1e-13 
                DegenerateStateIndices(NumDegenerateStates+1) = ll;
                NumDegenerateStates = NumDegenerateStates + 1;
            end
        end

        if NumDegenerateStates ~= kfArray(iFS,6)
            print('Number of degenerate FS does not match Number of degenerate Efs', NumDegenerateStates,kfArray(iFS,6));
        end
        
    for rr=1:NumDegenerateStates % look for the eigenvalue that matches the Fermi surface that we are on
        % todo: check that we should find only one eigenvalue that matches
        % teh Fermi surface.  Also,  maybe 
        % we could check whether the Fermi surface finder missed any solutions.

            r = DegenerateStateIndices(rr);

        % find this Fermi surface's contribution to the scalar density of
        % states.  This expression divides the Fermi momentum kf by vf.
        % Now, for the eigenvector whose eigenvalue is Ef, calculate the
        % expectation value of the velocity operator.  This is the Fermi
        % velocity at Ef and at angle thetak.
            % This calculates the velocity operator along the angle thetak.
            hhh2 = (hhha(:,:,2) * cos(thetak)) + (hhha(:,:,3) * sin(thetak));
        vvt = vv';
        thisvf = vvt(r,:) * hhh2 *  vv(:, r);

        DOSWeight = ((2 * pi/ NumValueTheta) * (kfArray( iFS,1) / thisvf ));
        DOSAtThetaArray(n) = DOSAtThetaArray(n) + DOSWeight;
    
            % calculate the matrix version of the density matrix
            RhoTmp = (vv(:,r) * vvt(r,:) ) * DOSWeight;
            RhoNArray(:,:,n) = RhoNArray(:,:,n) + RhoTmp;
            
            
            % calculate CollisionOperator
 CollisionOperator1 = zeros(LinearLength,LinearLength,Norbitals * Norbitals, Norbitals * Norbitals);
 DiffExpCoeffs1 = zeros(NumDiffExpCoeffs, Norbitals * Norbitals, Norbitals * Norbitals);
 
            % Four loops over all  matrix elements in the diffuson
            for s1 = 1:Norbitals
    for s2 = 1:Norbitals
        for s3 = 1:Norbitals
            for s4 = 1:Norbitals
                t1 = (s1 - 1) * 4 + s2;
                t2 = (s3 - 1) * 4 + s4;
                
                % calculate contributions to the matrix element from each
                % eigenstate.
                    for m = 1:Norbitals
        thisvfx = vvt(m,:) * hhha(:,:,2) *  vv(:, m) / hbarOverTau;
        thisvfy = vvt(m,:) * hhha(:,:,3) *  vv(:, m) / hbarOverTau;
%         thisvfx = 0;
%         thisvfy = 0;
                        delta_plus = dd(r,r) - dd(m,m);
                        mult1 =  vvt(r,s1) * vv(s3,r) * vv(s2,m) * vvt(m,s4) * DOSWeight * 0.5;
                        mult2 = conj(mult1); %todo
%   todo: figure this out        mult2 =  vvt(r,s4) * vv(s2,r) * vv(s3,m) * vvt(m,s1) * DOSWeight * 0.5;
                        ESOOverHbar = delta_plus / hbarOverTau;
                        a1 = 1 - (1i * ESOOverHbar);
                        a2 = 1 + (1i * ESOOverHbar);

% Coefficient for q-independent part of diffuson
DiffExpCoeffs1(1,t1,t2) = DiffExpCoeffs1(1,t1,t2) + (mult1 / a1) + (mult2/a2);
% Coefficient for (i qx)-dependent part of diffuson
DiffExpCoeffs1(2,t1,t2) = DiffExpCoeffs1(2,t1,t2) - (mult1 * thisvfx / (a1 * a1)) - (mult2 * thisvfx / (a2 * a2));
% Coefficient for (i qy)-dependent part of diffuson
DiffExpCoeffs1(3,t1,t2) = DiffExpCoeffs1(3,t1,t2) - (mult1 * thisvfy / (a1 * a1)) - (mult2 * thisvfy / (a2 * a2));
% Coefficient for (i qx)^2-dependent part of diffuson
DiffExpCoeffs1(4,t1,t2) = DiffExpCoeffs1(4,t1,t2) + (mult1 * thisvfx * thisvfx / (a1 * a1 * a1)) + (mult2 * thisvfx * thisvfx / (a2 * a2 * a2));
% Coefficient for (i qy)^2-dependent part of diffuson
DiffExpCoeffs1(5,t1,t2) = DiffExpCoeffs1(5,t1,t2) + (mult1 * thisvfy * thisvfy / (a1 * a1 * a1)) + (mult2 * thisvfy * thisvfy / (a2 * a2 * a2));
% Coefficient for (i qx)*(i qy)-dependent part of diffuson
DiffExpCoeffs1(6,t1,t2) = DiffExpCoeffs1(6,t1,t2) + (mult1 * thisvfx * thisvfy / (a1 * a1 * a1)) + (mult2 * thisvfx * thisvfy / (a2 * a2 * a2));

                        % a1= 1;
% a2 = 1;
                        % todo: implement cutoffs on ESOOverHbar(contribution is zero) and on
                        % KETauOverHbar (fixed contribution)
                        
        % loop over kspace grid
        for qx = 1:LinearLength 
                            for qy = 1:LinearLength
                                
%                                 qmag = kspacegrid(qx,qy,1);
%                                 costhetaq = kspacegrid(qx,qy,2);
%                                 sinthetaq = kspacegrid(qx,qy,3);
                                
        % Now, for the m-th eigenvector, calculate the
        % expectation value of the velocity operator.  This is the Fermi
        % velocity at the m-th Fermi surface and at angle thetak.
            % This next line calculates the velocity operator along the angle thetak.
            % this next line takes 7 % of the total computational time.
            KETauOverHbar =   kspacegrid(qx,qy,1) * ((thisvfx * kspacegrid(qx,qy,2)) + (thisvfy * kspacegrid(qx,qy,3))) ;
%             if(abs(real(KETauOverHbar)) > 0.000000000001)
%                 real(KETauOverHbar), imag(KETauOverHbar), qx, qy
%             end
            % This next line calculates the velocity along angle thetak.
            % This next line takes 67 % of the total computational time.
            delta1 =  (mult1 / (a1 + (1i *  KETauOverHbar))) + ( mult2 / (a2 + (1i * KETauOverHbar)));
 
                    % this next line takes 25 % of the total computational
                    % time.
                   CollisionOperator1(qx,qy,t1,t2) = CollisionOperator1(qx,qy,t1,t2) + delta1; 
                       
                            end % end of loops over qi, qy
                        end
                        
                    end % end of loop over orbitals

                    % end of four loops over all matrix elements in the diffuson
            end
        end 
    end
            end %last of the four loops
            
            % update the running total of CollisionOperator
             IntegratedCollisionOperator = IntegratedCollisionOperator + CollisionOperator1;
             DiffExpCoeffs = DiffExpCoeffs + DiffExpCoeffs1;
            
%        check whether multiplying by ProbVector gives RhoTmp
 res = zeros(Norbitals * Norbitals);
 resm = zeros(Norbitals, Norbitals);
% calculate the inverse
            res = squeeze(CollisionOperator1(1,1,:,:)) * ProbVector;            
            % res is a vector, now transfer it to resm which is a  matrix.
            for s1 = 1:Norbitals
                for s2 = 1:Norbitals
             t1 = (s2 - 1) * 4 + s1;
                    resm(s1, s2) = res(t1);
                end
            end
            % RhoTmp and resm should be absolutely identical because of the
            % way that RhoTmp, D_inverse1, and ProbVector were constructed.
            if norm(resm-RhoTmp) > 0.00000001
                %norm(resm-RhoTmp)
               % msgbox "probability not conserved."
            end
            % end of checking for probability conservation
            
            
    end  % end of loop over degenerate states
    
        
    end % end of loop/sum over Fermi surfaces



end % loop over theta

kfAtThetaZero = kfsArray(:,1); % For thetea = 0. Self note: Only pick first column since rest are repetition.
KfAtThetaPi2  = kfsArray(:,floor((NumValueTheta/ 4))); % For theta = approx pi/2.
% Integrate over theta
IntegratedDOS=sum(DOSAtThetaArray) ;
for n = 1: NumValueTheta  % loop over the angle around the origin
    IntegratedRhoN = IntegratedRhoN + RhoNArray(:,:,n);
end


% check whether the trace of IntegratedRhoN is equal to IntegratedDOS
if abs(trace(IntegratedRhoN) - IntegratedDOS) > 0.0000000001
    trace(IntegratedRhoN)
    IntegratedDOS
    msgbox 'the trace of IntegratedRhoN is not equal to IntegratedDOS';
end

%        check again whether multiplying by ProbVector gives Rho - see
%        comments on previous check
 res = zeros(Norbitals * Norbitals);
 resm = zeros(Norbitals, Norbitals);
            res = squeeze(IntegratedCollisionOperator(1,1,:,:)) * ProbVector;
            for s1 = 1:Norbitals
                for s2 = 1:Norbitals
             t1 = (s2 - 1) * 4 + s1;
                    resm(s1, s2) = res(t1);
                end
            end
            if norm(resm-IntegratedRhoN) > 0.00000001
                 norm(res), resm, IntegratedRhoN, resm-IntegratedRhoN
                msgbox '2 probability not conserved 2.'
            end
            %Normalize the Collision operator correctly.
IntegratedCollisionOperator = IntegratedCollisionOperator / (IntegratedDOS/Norbitals);
DiffExpCoeffs = DiffExpCoeffs / (IntegratedDOS/Norbitals);
% Dinv = 1 - I.
for qx = 1:LinearLength
    for qy = 1:LinearLength
IntegratedDinverse(qx,qy,:,:) = eye(Norbitals * Norbitals) - squeeze(IntegratedCollisionOperator(qx,qy,:,:));
    end
end
DiffExpCoeffs(1,:,:) = eye(Norbitals * Norbitals) - squeeze(DiffExpCoeffs(1,:,:));
DiffExpCoeffs(2:6,:,:) = - DiffExpCoeffs(2:6,:,:);

if norm(squeeze(IntegratedDinverse(1,1,:,:)) - squeeze(DiffExpCoeffs(1,:,:))) > 1e-10
    msgbox 'Zeroth order constant in diffusion expansion is wrong'
end
    
%  todo, handle these cases when LinearLength is even          if ((2 * (j-1)) == LinearLength) && ((2 * (i-1)) == LinearLength)
%               sinq = 0; cosq = 0;
%            elseif ((2 * (j-1)) == LinearLength) && i == 1
%               sinq = 0; cosq = 0;
%            elseif j == 1 && ((2 * (i-1)) == LinearLength)
%               sinq = 0; cosq = 0;
%            end
    
%squeeze(IntegratedDinverse(1,1,:,:))
% TestMatrix = zeros(LinearLength, LinearLength);
% for qx = 1:LinearLength
%     for qy = 1:LinearLength
%         TestMatrix(qx,qy) = ProbVector' * IntegratedDinverse(qx,qy,:,:) * ProbVector;
%     end
% end
% TestMatrix
% output
end