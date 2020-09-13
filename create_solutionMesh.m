%*********************************************************************
% Create the new mesh with including the new elements (points) and 
% ensures that the deleted elements do not appear anymore. 
% This code also proyect the last time step solution on the current ideal 
% mesh (required by the non-linear solver)
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Macrogeo,p_new] = create_solutionMesh(coord,element,...
    Nuevos_points,Delete_points,p_prev,sign)

global Lref
%% Macro_geo

Nuevos_points = unique(Nuevos_points,'rows');

if ~isempty(Delete_points)
Macrogeo.coordinate = setdiff(coord,Delete_points,'rows');
else
  Macrogeo.coordinate = coord;
end
Macrogeo.coordinate = unique([Macrogeo.coordinate;Nuevos_points],'rows');

Macrogeo.element    = delaunay(Macrogeo.coordinate(:,1),...
    Macrogeo.coordinate(:,2));
% nElement -> number of element at each mesh
Macrogeo.nElement  = size(Macrogeo.element,1);
% nnodes -> number of nodes at each mesh
Macrogeo.nnodes = size(Macrogeo.coordinate,1);

fprintf('\n New %i - Del %i \n',Macrogeo.nElement-size(element,1),size(Delete_points,1))

Macrogeo.size = [0 inf];
Macrogeo.bari = zeros(Macrogeo.nElement,2);

if sign ~= 0
    
    TR = triangulation(element,coord);
    p_new = zeros(Macrogeo.nElement,1);
    for j=1:Macrogeo.nElement
        
        % Coord (x;y) de cada uno de los vertices del tríangulo
        cc = Macrogeo.coordinate(Macrogeo.element(j,:),:)';
        
        Macrogeo.bari(j,:) = sum(cc')/3;
        elem = pointLocation(TR,Macrogeo.bari(j,:));
        
        bari = Macrogeo.bari(j,:);
    if isnan(elem)
        bari(1) = bari(1)+ 0.1*(1/Lref);
        elem = pointLocation(TR,bari);
    end
            p_new(j,1)= p_prev(elem);
    end
else
    p_new = zeros(Macrogeo.nElement,1);
    for j=1:Macrogeo.nElement
        cc = Macrogeo.coordinate(Macrogeo.element(j,:),:)';
        Macrogeo.bari(j,:) = sum(cc')/3;
    end
%     right = find(Macrogeo.bari(:,1)>0.5*(220/Lref));
%     p_new(Macrogeo.element(right,:),1)=1;
end

