%*********************************************************************
%*                                                                   *
%*                    Micro Boundary constraints                     *
%*                                                                   *
%*********************************************************************
%
% This function is usefull to impose the boundary conditions on the micro
% scale problem.
% Remark: The boundary condition will be Periodic or dirichlet, the user
% can handle it.
% Case Dirichlet: We impose weakle the Dirichlet boundary condition like an
% new data of the liner system.
% Periodic: We impose the peridicity of the soluton at the boundary using
% the selected pairs of edges (Pre-process) and we impose the condition
% like a constraints of the problem, for this reason we use 'lsqlin' to
% solve linear systems with constraints.
%
%***------------------------------------
%***Inputs: PosFlujo, edgesKnum (See Micro_solver)
%           H: Main matrix. J1/2: RHS of each cell problem.
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [SolFlujo1,SolFlujo2] = Micro_Boundary(M,b1,b2,Micro_geo)

nEdgeTotal   = Micro_geo.noedges;
nElement     = Micro_geo.nElement;
edgePar      = Micro_geo.edgePar;

IndicadorDirichFlujo = eye(nEdgeTotal+nElement);
for kk = 1:size(edgePar)
    IndicadorDirichFlujo(edgePar(kk,2),edgePar(kk,1)) = -1;
    IndicadorDirichFlujo(nEdgeTotal+edgePar(kk,4),nEdgeTotal+edgePar(kk,3)) = 1;
end
IndicadorDirichFlujo(:,[edgePar(:,2);nEdgeTotal+edgePar(:,4)])=[];
% IndicadorDirichFlujo(:,nEdgeTotal+edgePar(:,4))=[];

MM   = (IndicadorDirichFlujo'*M*IndicadorDirichFlujo);
b1 = (IndicadorDirichFlujo'*b1);
b2 = (IndicadorDirichFlujo'*b2);

posDir  = [1 2]; %edge 1 imposed to be 0
freePOS = setdiff(1:length(b1),posDir);

% -----------------------------------------------------------
%            SOLUTION OF DE EDGE PROBLEM
% -----------------------------------------------------------

Sol1(posDir,1)  = 0;
Sol1(freePOS,1) = MM(freePOS,freePOS)\b1(freePOS);
%         Sol1(freePOS,1) = gmres(HH(freePOS,freePOS),JJ11(freePOS));
Sol2(posDir,1)  = 0;
Sol2(freePOS,1) = MM(freePOS,freePOS)\b2(freePOS);

% Copy the solution
indSolut = setdiff(1:nEdgeTotal+nElement,[edgePar(:,2);nEdgeTotal+edgePar(:,4)]);
% indSolut = setdiff(1:nEdgeTotal+nElement,nEdgeTotal+edgePar(:,4));
SolFlujo1 = zeros(nEdgeTotal+nElement,1);
SolFlujo2 = zeros(nEdgeTotal+nElement,1);

SolFlujo1(indSolut,1)    = Sol1;
% SolFlujo1([edgePar(:,2);nEdgeTotal+edgePar(:,4)],1) = ...
%     SolFlujo1([edgePar(:,1);nEdgeTotal+edgePar(:,3)],1);
SolFlujo1(edgePar(:,2),1) =  -SolFlujo1(edgePar(:,1),1);
SolFlujo1(nEdgeTotal+edgePar(:,4),1) =  SolFlujo1(nEdgeTotal+edgePar(:,3),1);

SolFlujo2(indSolut,1)    = Sol2;
SolFlujo2(edgePar(:,2),1) =  -SolFlujo2(edgePar(:,1),1);
SolFlujo2(nEdgeTotal+edgePar(:,4),1) =  SolFlujo2(nEdgeTotal+edgePar(:,3),1);






