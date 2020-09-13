%*********************************************************************
%*                                                                   *
%*              EFECTIVE Diffusion tensor K _ Upscaled               *
%*                                                                   *
%*********************************************************************
%
% This functions compute the efective difussion tensor if the original
% tensor does not dependent of the macro scale.
%
%***------------------------------------
%***Inputs: Solution of the micro cell problems.
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [K_Efective] = EfectivePermTensor(Micro_geo,xx,yy,K_micro11,K_micro12,K_micro21,K_micro22,Vel1,Vel2)

element      = Micro_geo.element;
coordinate   = Micro_geo.coordinate;

AreaSum    = 0;
K_Efective = zeros(2,2);

for j = 1:Micro_geo.nElement
    % Coord (x;y) triangle vertex
    coord = coordinate(element(j,:),:)';
    
    % coord egdes : [1:2] inicial , [3:4] final
    p      = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    % Length edge
    le     = [norm(p(1:2,1)-p(3:4,1)) norm(p(1:2,2)-p(3:4,2))...
        norm(p(1:2,3)-p(3:4,3))];
    
    I = diag(Micro_geo.nodes2edge(Micro_geo.element(j,[2 3 1]),...
        Micro_geo.element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==Micro_geo.edge2element(I,4)))= -1;
    
    bari = sum(coord,2)/3;
    area = det([1,1,1;coord])/2;
    
    N=bari(:)*ones(1,3)-coord;
    P1=N*1/2*diag(signum)*diag(le)*Vel1(I);
    P2=N*1/2*diag(signum)*diag(le)*Vel2(I);
    
    aux11 = interp2(xx,yy,K_micro11',bari(1),bari(2),'nearest');
    aux12 = interp2(xx,yy,K_micro12',bari(1),bari(2),'nearest');
    aux21 = interp2(xx,yy,K_micro21',bari(1),bari(2),'nearest');
    aux22 = interp2(xx,yy,K_micro22',bari(1),bari(2),'nearest');
    
    K_aux = [aux11 aux12; aux21 aux22];
    
    K_Efective = K_Efective + (K_aux*area +[P1 P2]);
    
    % A_EFECTIVE = INT_P (e_j+w^j)*e_i
    AreaSum = AreaSum + area;  
end
K_Efective =  K_Efective./AreaSum;