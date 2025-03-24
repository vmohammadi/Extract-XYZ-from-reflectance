function [deltaElabD65] = XYZ2DeltaElab(XYZ_r,XYZ_Pr,illuminant,observer)
  %% This function calculates the delta E LAB from two sets of XYZ values

  %            XYZ_r  :   xyz from Actual Reflectance that is calculated by known illuminat and observer (3,1)  
  %           XYZ_Pr  :   predicted xyz under known illuminat and observer (3,1)  
  %       illuminant  :   SPD of illuminant        (n*2)     ,   first coloumn  belongs to  the wavelengths(nm)
  %         observer  :   color matching functions (n*4)     ,   first coloumn  belongs to  the wavelengths(nm)
 
    interval    =        10     %  (nm)       interval  :   interval between wavelength
%%  k is the regularization factor 
% load  color matching functions of 1931 for 2 degree observer

    load ('ciedeg2.txt');

    observer   =     ciedeg2(41:interval:341,:); % cut the part of spectrum between 400 and 700
    illuminant = illuminant(9:2:69,:);
      k   =     100/sum(observer(:,3).* illuminant(:,2));

%%   XYZn is XYZ of illuminant

    [m,n] = size(observer);
    B = zeros(m,n);
    B(:,1) = observer(:,1);
    
      for  i=1:3
          
          B(:,i+1)   =    k*observer(:,i+1).*illuminant(:,2);
      end  
      
      % for Sample
             Xn      =   sum(sum(bsxfun(@times,B(:,2),XYZ_r(:,1)'))) ; 
             Yn      =   sum(sum(bsxfun(@times,B(:,3),XYZ_r(:,1)'))) ; 
             Zn      =   sum(sum(bsxfun(@times,B(:,4),XYZ_r(:,1)'))) ; 

       XYZn    =    [sum(Xn);sum(Yn);sum(Zn)];    %  size = (3*1)

%% lab_r from XYZ_r

        for i=1:3
             
              if (XYZ_r(1,:)/Xn)>0.008856
                     L_r   =    (116*(XYZ_r(2,:)/Yn)^(1/3))-16  ;
                     a_r   =    500*((XYZ_r(1,:)/Xn)^(1/3)-(XYZ_r(2,:)/Yn)^(1/3)) ; 
                     b_r   =    200*((XYZ_r(2,:)/Yn)^(1/3)-(XYZ_r(3,:)/Zn)^(1/3)) ;
              elseif XYZ_r(2,:)/Yn>0.008856
                     L_r   =    (116*(XYZ_r(2,:)/Yn).^(1/3))-16  ;
                     a_r   =    500*((XYZ_r(1,:)/Xn).^(1/3)-(XYZ_r(2,:)/Yn).^(1/3)) ; 
                     b_r   =    200*((XYZ_r(2,:)/Yn).^(1/3)-(XYZ_r(3,:)/Zn).^(1/3)) ;
              elseif XYZ_r(3,:)/Zn>0.008856
                     L_r   =    (116*(XYZ_r(2,:)/Yn).^(1/3))-16  ;
                     a_r   =    500*((XYZ_r(1,:)/Xn).^(1/3)-(XYZ_r(2,:)/Yn).^(1/3)) ; 
                     b_r   =    200*((XYZ_r(2,:)/Yn).^(1/3)-(XYZ_r(3,:)/Zn).^(1/3)) ;
              else
                     L_r   =    116*7.787*(XYZ_r(2,:)/Yn);
                     a_r   =    500*7.787*((XYZ_r(1,:)/Xn)-(XYZ_r(2,:))/Yn) ; 
                     b_r   =    200*7.787*((XYZ_r(2,:)/Yn)-(XYZ_r(3,:)/Zn)) ;
             end
        end

                     Lab_r    =    [L_r;a_r;b_r]

%% lab_Pr from XYZ_Pr

        for i=1:3
             
              if XYZ_Pr(1,:)/Xn>0.008856
                     L_Pr   =    (116*(XYZ_Pr(2,:)/Yn).^(1/3))-16  ;
                     a_Pr   =    500*((XYZ_Pr(1,:)/Xn).^(1/3)-(XYZ_Pr(2,:)/Yn).^(1/3)) ; 
                     b_Pr   =    200*((XYZ_Pr(2,:)/Yn).^(1/3)-(XYZ_Pr(3,:)/Zn).^(1/3)) ;
              elseif XYZ_Pr(2,:)/Yn>0.008856
                     L_Pr   =    (116*(XYZ_Pr(2,:)/Yn).^(1/3))-16  ;
                     a_Pr   =    500*((XYZ_Pr(1,:)/Xn).^(1/3)-(XYZ_Pr(2,:)/Yn).^(1/3)) ; 
                     b_Pr   =    200*((XYZ_Pr(2,:)/Yn).^(1/3)-(XYZ_Pr(3,:)/Zn).^(1/3)) ;
              elseif XYZ_Pr(3,:)/Zn>0.008856
                     L_Pr   =    (116*(XYZ_Pr(2,:)/Yn).^(1/3))-16  ;
                     a_Pr   =    500*((XYZ_Pr(1,:)/Xn).^(1/3)-(XYZ_Pr(2,:)/Yn).^(1/3)) ; 
                     b_Pr   =    200*((XYZ_Pr(2,:)/Yn).^(1/3)-(XYZ_Pr(3,:)/Zn).^(1/3)) ;              
              else
                     L_Pr   =    116*7.787*(XYZ_Pr(2,:)/Yn);
                     a_Pr   =    500*7.787*((XYZ_Pr(1,:)/Xn)-(XYZ_Pr(2,:))/Yn) ; 
                     b_Pr   =    200*7.787*((XYZ_Pr(2,:)/Yn)-(XYZ_Pr(3,:)/Zn)) ;
             end
        end

                     Lab_Pr    =    [L_Pr;a_Pr;b_Pr]
         
%% Delta_E calculation between Lab_r & LabPr

                DeltaE1976  =   sqrt(sum((Lab_r-Lab_Pr).^2))
                deltaElabD65 = DeltaE1976;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
