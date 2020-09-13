
%*********************************************************************
% Aporteriori refinement
% This code use Aposteriori.m 
%*********************************************************************
%
%***------------------------------------
%***Inputs: Heterogeneos permeability field
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [A_refPermeability,New_points,Delete_points,Ind_refin]=...
    AposterioriRefinement(eta_T,Klevels,A_refPermeability,Ind_refin,field,time)

global Macro_geo N_real Num_levels Lref
global crit_ref crit_coar

New_points = [];
Delete_points = [];

% Refinement criteria

x_ref = 0:N_real(1); y_ref = 0:N_real(2);
x_ref = x_ref/Lref;  y_ref = y_ref/Lref;

xbin = discretize (Macro_geo.(field).bari(:,1),x_ref);  % de que columna?
ybin = discretize (Macro_geo.(field).bari(:,2),y_ref);  % de que fila?
idx = sub2ind(size(Ind_refin),xbin,ybin);
levels = Ind_refin(idx);

criterio   = crit_ref*max(eta_T);
criterio0  = crit_coar*min(eta_T);
[REFINAR] = find(eta_T > criterio);
[COARSEN] = find(eta_T < criterio0);

eval_level= Macro_geo.(field).bari(REFINAR,:);

% Auxiliar matrix to no refine twice
Aux_Ind_refin = zeros(size(Ind_refin));

%%
for j = 1:length(REFINAR)
    xbin = discretize (eval_level(j,1),x_ref);  % de que columna?
    ybin = discretize (eval_level(j,2),y_ref); % de que fila?
    
    level = Ind_refin(xbin,ybin);
    
    if ~isnan(Aux_Ind_refin(xbin,ybin))
        
        if level ~= Num_levels
            
            f_level  = sprintf('ref%i',level);
            f_level1 = sprintf('ref%i',level+1);
            
            mesh_level = Klevels.(f_level).gridmesh{1};
            xxedge = linspace(0,N_real(1)/Lref,size(mesh_level,2));
            yyedge = linspace(0,N_real(2)/Lref,size(mesh_level,1));
            xbin_level = discretize (eval_level(j,1),xxedge);
            ybin_level = discretize (eval_level(j,2),yyedge);
            
            % NewPoints
            
            new12= [ones(2,1)*(xxedge(xbin_level) + 1/2*(xxedge(xbin_level+1)-xxedge(xbin_level))),...
                [yyedge(ybin_level);yyedge(ybin_level+1)]];
            new34=[ [xxedge(xbin_level);xxedge(xbin_level+1)],...
                ones(2,1)*(yyedge(ybin_level) + 1/2*(yyedge(ybin_level+1)-yyedge(ybin_level)))];
            new5 = [ (xxedge(xbin_level) + 1/2*(xxedge(xbin_level+1)-xxedge(xbin_level))),...
                (yyedge(ybin_level) + 1/2*(yyedge(ybin_level+1)-yyedge(ybin_level)))];
            New_points = [New_points; new12; new34; new5];
            
            %%
            % Duplicar valores, malla refinada (artificial)
            refinar_c = [2*xbin_level-1; 2*xbin_level-1; 2*xbin_level; 2*xbin_level];
            refinar_r = [2*ybin_level-1; 2*ybin_level; 2*ybin_level-1; 2*ybin_level];
            
            mesh_level1 = Klevels.(f_level1).gridmesh{1};
            ref_index = unique(sub2ind(size(mesh_level1),refinar_r,refinar_c));
            
            % malla refinada (artificial)
            aux1= Klevels.(f_level1).tensor(:,:,1)';
            aux2= Klevels.(f_level1).tensor(:,:,2)';
            aux3= Klevels.(f_level1).anisot(:,:,1)';
            aux4= Klevels.(f_level1).anisot(:,:,2)';
            valores1 = zeros(2);
            valores2 = zeros(2);
            valores3 = zeros(2);
            valores4 = zeros(2);
            valores1(:) = aux1(ref_index); %valores
            valores2(:) = aux2(ref_index);
            valores3(:) = aux3(ref_index);
            valores4(:) = aux4(ref_index);
            
            if level+1 <Num_levels
                valores1 = kron(valores1,ones(2*(Num_levels-level-1)));
                valores2 = kron(valores2,ones(2*(Num_levels-level-1)));
                valores3 = kron(valores3,ones(2*(Num_levels-level-1)));
                valores4 = kron(valores4,ones(2*(Num_levels-level-1)));
            end
            
            aay = 2*(Num_levels-level)*(xbin_level-1)+1:2*(Num_levels-level)*xbin_level;
            aax = 2*(Num_levels-level)*(ybin_level-1)+1:2*(Num_levels-level)*ybin_level;
            %             ref_index = unique(sub2ind([N_real 2],aay,aax,ones(size(aay))));
            
            A_refPermeability(aay,aax,1) = valores1';
            A_refPermeability(aay,aax,2) = valores2';
            A_refPermeability(aay,aax,3) = valores3';
            A_refPermeability(aay,aax,4) = valores4';
            
            Ind_refin(aay,aax) = level+1;
            %             Si_Refinar = [Si_Refinar;REFINAR(j)];
            Aux_Ind_refin(aay,aax) = NaN;
            
        end
    end
end

%%
eval_level= Macro_geo.(field).bari(COARSEN,:);
for j = 1:length(COARSEN)
    xbin = discretize (eval_level(j,1),x_ref);  % de que columna?
    ybin = discretize (eval_level(j,2),y_ref); % de que fila?
    
    level = Ind_refin(xbin,ybin);
    
    if ~isnan(Aux_Ind_refin(xbin,ybin))
        
        if level > 0
            
            f_level1  = sprintf('ref%i',level);
            f_level = sprintf('ref%i',level-1);
            
            %Grid de la malla que quiero
            mesh_level = Klevels.(f_level).gridmesh{1};
            xxedge = linspace(0,N_real(1)/Lref,size(mesh_level,2));
            yyedge = linspace(0,N_real(2)/Lref,size(mesh_level,1));
            xbin_level = discretize (eval_level(j,1),xxedge);
            ybin_level = discretize (eval_level(j,2),yyedge);
            
            % Delete points
            new12= [ones(2,1)*(xxedge(xbin_level) + 1/2*(xxedge(xbin_level+1)-xxedge(xbin_level))),...
                [yyedge(ybin_level);yyedge(ybin_level+1)]];
            new34=[ [xxedge(xbin_level);xxedge(xbin_level+1)],...
                ones(2,1)*(yyedge(ybin_level) + 1/2*(yyedge(ybin_level+1)-yyedge(ybin_level)))];
            new5 = [ (xxedge(xbin_level) + 1/2*(xxedge(xbin_level+1)-xxedge(xbin_level))),...
                (yyedge(ybin_level) + 1/2*(yyedge(ybin_level+1)-yyedge(ybin_level)))];
            Delete_points = [Delete_points;new12; new34; new5];
            
            %%
            % Duplicar valores, malla refinada (artificial)
%             refinar_c = [2*xbin_level-1; 2*xbin_level-1; 2*xbin_level; 2*xbin_level];
%             refinar_r = [2*ybin_level-1; 2*ybin_level; 2*ybin_level-1; 2*ybin_level];
            
%             mesh_level1 = Klevels.(f_level1).gridmesh{1};
%             ref_index = unique(sub2ind(size(mesh_level1),refinar_r,refinar_c));
            
            % malla refinada (artificial)
            aux1= Klevels.(f_level).tensor(:,:,1)';
            aux2= Klevels.(f_level).tensor(:,:,2)';
            aux3= Klevels.(f_level).anisot(:,:,1)';
            aux4= Klevels.(f_level).anisot(:,:,2)';
            valores1 = zeros(2);
            valores2 = zeros(2);
            valores3 = zeros(2);
            valores4 = zeros(2);
            valores1(:) = aux1(ybin_level,xbin_level); %valores
            valores2(:) = aux2(ybin_level,xbin_level);
            valores3(:) = aux3(ybin_level,xbin_level);
            valores4(:) = aux4(ybin_level,xbin_level);
            
            
            if level <Num_levels
                valores1 = kron(valores1,ones(2*(Num_levels-level)));
                valores2 = kron(valores2,ones(2*(Num_levels-level)));
                valores3 = kron(valores3,ones(2*(Num_levels-level)));
                valores4 = kron(valores4,ones(2*(Num_levels-level)));
            end
            
            aay = 2*(Num_levels-level+1)*(xbin_level-1)+1:2*(Num_levels-level+1)*xbin_level;
            aax = 2*(Num_levels-level+1)*(ybin_level-1)+1:2*(Num_levels-level+1)*ybin_level;
            %             ref_index = unique(sub2ind([N_real 2],aay,aax,ones(size(aay))));
            
            %             aay = xbin_level:xbin_level+1;
            %             aax = ybin_level:ybin_level+1;
            
            A_refPermeability(aay,aax,1) = valores1;
            A_refPermeability(aay,aax,2) = valores2;
            A_refPermeability(aay,aax,3) = valores3;
            A_refPermeability(aay,aax,4) = valores4;
            
            Ind_refin(aay,aax) = level-1;
            %             Si_Refinar = [Si_Refinar;REFINAR(j)];
            Aux_Ind_refin(aay,aax) = NaN;
            
        end
    end
end


% New_points = unique(New_points,'rows');
