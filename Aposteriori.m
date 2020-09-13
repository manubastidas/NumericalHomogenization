%*********************************************************************
% Aposteriori error estimation
%
% This code is based on:
% Bahriawati, C., & Carstensen, C. (2005). 
% Three MATLAB implementations of the lowest-order Raviart-Thomas 
% MFEM with a posteriori error control. 
% Computational Methods in Applied Mathematics, 5(4), 333-361.
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium


function eta_T = Aposteriori(element,coordinate,pEval)
             
u_h = zeros(size(coordinate,1),2);
supp_area = zeros(size(coordinate,1),1);

% SURRONDING AREA
for j = 1 : size(element,1)
   supp_area(element(j,:)) = supp_area(element(j,:)) + ...
   ones(3,1)*det([1,1,1;coordinate(element(j,:),:)'])/6;
   u_h(element(j,:),:) = u_h(element(j,:),:) + ...
         det([1,1,1;coordinate(element(j,:),:)']) * ...
         ( (pEval(3*(j-1)+[1,2,3],:))' * [4 1 1;1 4 1;1 1 4]/36)';
end
% AVERAGING AREA 
u_h=u_h./(supp_area*ones(1,2));

% L2 NORM
eta_T=zeros(size(element,1),1);
for j = 1 : size(element,1)  
  eta_T(j) = sqrt(det([1,1,1;coordinate(element(j,:),:)']) * ...
              (sum( ( [4 1 1;1 4 1;1 1 4]/6 * u_h(element(j,:),:) -...
               pEval(3*(j-1)+[1,2,3],:)).^2') * ones(3,1)/6 ));           
end        
% eta_T = sqrt(sum(eta_T.^2));  
 

















