%*********************************************************************
%*                        L2 ERROR MFEM
%*********************************************************************
%
% This function calculate the L2 error of the FEM aproximation.
% L2 error --> SQRT (sum(K in trian) int_K (aprox - exact)^2)
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Error_pTCont, Error_p, abs_diff] = errorL2_Macro(p,pExacta,geo,ref_geo)

global Lref
% structure for the coarse mesh
TR = triangulation(geo.element,geo.coordinate);

Error_pT     = zeros(size(ref_geo.element,1),1);
Error_pTCont = zeros(size(ref_geo.coordinate,1),1);
abs_diff     = zeros(size(ref_geo.coordinate,1),1);

contador = zeros(size(ref_geo.coordinate,1),1);

for j=1:ref_geo.nElement
    
    % Coord (x;y) de cada uno de los vertices del tríangulo
    coord = ref_geo.coordinate(ref_geo.element(j,:),:)';
    
    % baricenter of the coarse element
    bari = sum(coord,2)/3;
    coarse_elem = pointLocation(TR,bari');
    
    if isnan(coarse_elem)
        bari(1) = bari(1)+ 0.1*(1/Lref);
        coarse_elem = pointLocation(TR,bari');
    end
    diff = (p(3*coarse_elem)-pExacta(3*j));
    
    abs_diff(ref_geo.element(j,:)) = abs_diff(ref_geo.element(j,:)) + abs(diff);
    Error_pT(j) = diff.^2*det([1,1,1;coord])/2;
    Error_pTCont(ref_geo.element(j,:)) = Error_pTCont(ref_geo.element(j,:))+ Error_pT(j);
    
    contador(ref_geo.element(j,:)) = contador(ref_geo.element(j,:)) +1;
end
abs_diff     = abs_diff./contador;
Error_pTCont = Error_pTCont./contador;

% abs_diff     = max(eps,abs_diff)./contador;
% Error_pTCont = max(eps,Error_pTCont)./contador;
Error_p  = sqrt(sum(Error_pT));
