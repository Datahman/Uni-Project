    function [excitonkspace1] = excitonkspace(LatticeSpacingOverUnitLength, LinearLength, ThetaX)
    
    excitonkspace1 = zeros(LinearLength, LinearLength, 3);
   % qi,qj go  up to 2 pi / LatticeSpacingOverUnitLength in increments of 2 pi / (LinearLength * LatticeSpacingOverUnitLength)   
   % i=1 is the same as zero momentum
   for i=1:LinearLength
       for j=1:LinearLength
                    
               %the next two lines make the dispersion periodic on the
               %Brillioun zone
               qi = sin((i-1) *  pi / LinearLength) * (2 / LatticeSpacingOverUnitLength);
               qj = sin((j-1) *  pi / LinearLength) * (2 / LatticeSpacingOverUnitLength);
               qlsqr = qi^2 + qj^2; 

              % The following computes sin(N(thetaq- ThetaX)) and
               % cos(N(thetaq-ThetaX))
               costhetaq = 0; sinthetaq = 0;
            if i == 1 && j == 1 
               sinq = 0; cosq = 0; % not used because we already set costhetaq=sinthetaq=0
           elseif ((2 * (j-1)) == LinearLength) && ((2 * (i-1)) == LinearLength)
              sinq = 0; cosq = 0;
           elseif ((2 * (j-1)) == LinearLength) && i == 1
              sinq = 0; cosq = 0;
           elseif j == 1 && ((2 * (i-1)) == LinearLength)
              sinq = 0; cosq = 0;
            else
                qix = sin((i-1) *  2 *  pi / LinearLength) * (1 / LatticeSpacingOverUnitLength);
               qjx = sin((j-1) *  2 *  pi / LinearLength) * (1 / LatticeSpacingOverUnitLength);
               qlsqrx = qix^2 + qjx^2; 
               sinq = qix / sqrt(qlsqrx);
               cosq = qjx / sqrt(qlsqrx);
            
              costhetaq = ((cosq * cos(ThetaX)) + (sinq * sin(ThetaX)) ) ;
              sinthetaq = ((-cosq * sin(ThetaX)) + (sinq * cos(ThetaX)) ) ;
            end
    % todo: maybe fix the problem when qlsqr=0 and then we take its square root
    % however matlab is not throwing an exception, and the sqrt(0) is
    % multiplied by sinthetaq and costhetaq which are also zero.
               excitonkspace1(i,j,:) = [sqrt(qlsqr), sinthetaq, costhetaq];

       end % end of loop over j, which is the y coordinate
   end % end of loop over i, which is the x coordinate
   
    end %function
     
