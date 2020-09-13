
%*********************************************************************
% Main Code for Numerical adaptive homogenization
% Non-periodic case
%*********************************************************************
%
%***------------------------------------
%***Inputs: Heterogeneos permeability field
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

clear
close all
clc

%% Load the SPE10 DATA
name = 'SPEtop';
compute_Ref = 1;
compute_lev = 1;

addpath('GraphicsSPE10')
% file    = 'kk_b';
% field1  = load(file);
% field   = field1(2:2:end,1);
% K_perm  = reshape(field,60,220);

load('spe10_rock.mat')
file    = rock.perm(:,1);
K_perm_All  = reshape(file,60,220,85);
K_perm = K_perm_All(:,:,1);

%% INPUTS
global  N_real Macro_geo Time Num_levels
global b_nolin L_scheme bp_nolin nondimC
global Lref Harm ISO
global crit_ref crit_coar

N_real  = [220,60];
% N_real  = [resolution,resolution];
Harm = 0; %1 MEANS HARMONIC
ISO  = 1;

T_MAX = 1;
Time.tnSteps      = 10;
Time.dt           = T_MAX/Time.tnSteps;
Time.time_vec     = 0+Time.dt:Time.dt:T_MAX;

K_perm  = K_perm(1:N_real(2),1:N_real(1));
K_perm  = [K_perm K_perm(:,end)];
K_perm  = [K_perm;K_perm(end,:)];
% K_perm = K_perm./(max(max(K_perm)));

K_perm = K_perm*9.86E-16;
maxKperm = max(max(K_perm));

Lref = 33;
Tref = 1E3;
Pref = 1;
Cte  = 1;
alpha = (Lref^2*Pref^2)/(Tref*maxKperm);
nondimC = 1E-4;

K_perm = K_perm./maxKperm;

% Computing Isotropic vs Anisotropic
if ISO == 0
    %     rot = [cosd(25) -sind(25); sind(25) cosd(25)];
    %     Mmat = [1 0; 0 10^-3];
    rot = eye(2);
    Mmat = [1 0; 0 1];
    K_perm_full = zeros(N_real(2)+1,N_real(1)+1,4);
    for ri = 1:N_real(2)+1
        for rj = 1:N_real(1)+1
            K_aux = rot*K_perm(ri,rj)*Mmat*inv(rot);
            K_perm_full(ri,rj,1) = K_aux(1,1);
            K_perm_full(ri,rj,2) = K_aux(2,2);
            K_perm_full(ri,rj,3) = K_aux(1,2);
            K_perm_full(ri,rj,4) = K_aux(2,1);
        end
    end
    K_perm = K_perm_full;
else
    K_perm_full(:,:,1) = K_perm;
    K_perm_full(:,:,2) = K_perm;
    K_perm_full(:,:,3) = 0*K_perm;
    K_perm_full(:,:,4) = 0*K_perm;
    K_perm = K_perm_full;
end

NCoarse = [55 15]; % SPE10

% Parameters for the refinement
crit_ref = 0.5;
crit_coar = 10;

CurrentDir = pwd();
resultDir = [CurrentDir,'\',name];
mkdir(resultDir)

%% Reference Solution

save_file0 = [resultDir,'\','Reference_prev_',name,'.mat'];
if compute_Ref == 1
    [Macro_geoREF,Macro_SolREF,A,ele1,ele2,T] = ReferenceSolution_FEM_PREV(K_perm);
    save(save_file0,'Macro_geoREF','Macro_SolREF','A','ele1','ele2','T');
else
    load(save_file0)
    compute_Ref = 0;
end

b_nolin  = @(p) (p).^3*nondimC;
L_scheme = 3*nondimC;
bp_nolin = @(p) 3.*p.^2*nondimC;

save_file1 = [resultDir,'\','Reference_',name,'.mat'];
if compute_Ref == 1
    [Macro_SolREF] = ReferenceSolution_FEM(Macro_geoREF,Macro_SolREF,A,ele1,ele2,T);
    save(save_file1,'Macro_SolREF')
else
    load(save_file1)
    compute_Ref = 0;
end

%% Permeability levels
% CREAR LOS NIVELES DE LA PERMEABILIDAD
% upscaled permeabilities = A_Efective
% numl = numero de niveles que hemos creado

tic;
Nmicro  = N_real./NCoarse;

save_file2 = [resultDir,'\','Permlevels_',name,'.mat'];
if compute_lev == 1
    [Klevels,Num_levels]   = permeability_levels(K_perm,NCoarse);
    [Klevels_EfectiveFull] = permeability_levels2(Klevels,Num_levels);
    time_perm = toc;
    save(save_file2,'Klevels','Num_levels','Klevels_EfectiveFull','time_perm')
else
    load(save_file2)
    compute_lev = 0;
end


%% Solution + Aposteriori + Refinement
% REFINAR LA PERMEABILIDAD
% -----------------------------------------------------
% A_EfectivePerm = permeabilidad refinada en cada nivel

Macro_geo      = struct();
Macro_Sol      = struct();
EfectivePerm_Refin = struct();

Error_pT  = cell(0);
Error_p   = ones(Time.tnSteps,1);
Pressdiff = cell(0);

[xx,yy] = meshgrid(0:Nmicro(1):N_real(1),0:Nmicro(2):N_real(2));
xx = xx/Lref; yy=yy/Lref;
field = sprintf('time%i',1);

% -------------------------TIME 0-------------------------
% Solution 0 - Coarse mesh - KEffective Coarse
[Macro_geo.(field),p_prev] = create_solutionMesh([xx(:),yy(:)],[],[],[],[],0);

EfectivePerm_Refin.(field) = Klevels_EfectiveFull.ref0;

% Indices de refinamiento
Ind_refin  = zeros(N_real(1)+1,N_real(2)+1,Time.tnSteps+2);

Aux_points = [];

tic
% ------------------------TIME 1 to T-----------------------
for kk = 1:Time.tnSteps
    
    fprintf('\n Time %i/%i',kk,Time.tnSteps);
    
    %% MFEM Solution
    [A,elegidos,T]  = FEM_solutionPreprocess(field);
    
    [Macro_Sol.(field),p_prev] = FEM_solution(A,T,elegidos{1},elegidos{2},...
        EfectivePerm_Refin,field,p_prev,kk);
    
    Macro_Sol.(field).Mag_grad = sqrt(sum(Macro_Sol.(field).VelCont,2).^2);
    
    %%  NEXT PERMEABILITY & MESH
    field1 = sprintf('time%i',kk+1);
    
    if crit_ref ~= 0
        % Aposteriori error estimator
        eta_T = Aposteriori(Macro_geo.(field).element,...
            Macro_geo.(field).coordinate,Macro_Sol.(field).Vel);
        
        % --------------- Refine mesh
        [EfectivePerm_Refin.(field1),New_points,Delete_points,Ind_refin(:,:,kk+1)] = ...
            AposterioriRefinement(eta_T,Klevels,...
            EfectivePerm_Refin.(field),Ind_refin(:,:,kk),field,kk);
        
        [Macro_geo.(field1),p_prev] = create_solutionMesh...
            (Macro_geo.(field).coordinate,Macro_geo.(field).element,...
            New_points,Delete_points,p_prev,1);
    else
        Macro_geo.(field1) = Macro_geo.(field);
        EfectivePerm_Refin.(field1) = EfectivePerm_Refin.(field);
    end
    
    %% ERRORS
    [Error_pT{kk}, Error_p(kk), Pressdiff{kk}]  =...
        errorL2_Macro(Macro_Sol.(field).Pres,...
        Macro_SolREF.(field).Pres,Macro_geo.(field),Macro_geoREF);
    vElement(kk) = Macro_geo.(field).nElement;
    %
    field = field1;
end

%%
Error_T = sqrt(Time.dt*sum((Error_p.^2)));

% SAVING THE RESULTS
format shorte
fprintf('\n Error %i \n',Error_T)
format short
fprintf(' Av.Elements %i \n',ceil(mean(vElement)))
fprintf(' Last Time %i \n',vElement(end))


trin = sprintf('Solution_cr%i_cc%i_',crit_ref*10,crit_coar);
save_file = [resultDir,'\',trin,name,'.mat'];
time_solution = toc;

save(save_file)

