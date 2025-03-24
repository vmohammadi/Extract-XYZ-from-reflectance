function [XYZ] = R2XYZ(ActualR,illuminant,observer )
%%  This function extracts the XYZ values from Reflectance

  %         ActualR  :    Actual reflectance      (n*2)     , first coloumn  belongs to  the wavelengths(nm)
  %       illuminant :   SPD of illuminant        (n*2)     , first coloumn  belongs to  the wavelengths(nm)
  %         observer :   color matching functions (n*4)     , first coloumn  belongs to  the wavelengths(nm)
  
% All wavelengths  must be the same
% All date must be divide into the data at 560 nm to achieve the relative data


%%  k regularization factor 
% load  color matching functions of 1931 2 degree

      load ('ciedeg2.txt');

    observer   =     ciedeg2(41:10:341,:);
    
           k   =     1/sum(observer(:,3).* illuminant(:,2));

%%   XYZ calculation
    [m,n] = size(observer);
    B = zeros(m,n);
    B(:,1) = observer(:,1);
    
      for  i=1:3
          
          B(:,i+1)=k*observer(:,i+1).*illuminant(:,2);
      end  
      
      % for Sample
       X1    =   bsxfun(@times,B(:,2),ActualR(:,1)) ;  
       Y1    =   bsxfun(@times,B(:,3),ActualR(:,1)) ;  
       Z1    =   bsxfun(@times,B(:,4),ActualR(:,1)) ; 

       
       XYZ    =    [sum(X1);sum(Y1);sum(Z1)];    %  size = (3*1)

   
