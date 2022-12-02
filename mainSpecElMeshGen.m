%---------------------------------------------------------------
% NURBS-Enhanced Mesh Generator for
% Spectral Element Method
%---------------------------------------------------------------
%
% Read the Geometry imported from Rhino:
FileName = 'specElMesh_sample1_';
numPatch = 1;
%Degrees of Freedom per each node:
local_dof = 1;
%CREATE 2D IGA MESH (reads FileName):
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
%---------------------------------------------------------------
% Plot Imported 2-D NURBS Structure
iga2DmeshPlotNURBS(Nurbs2D);
%---------------------------------------------------------------
% Sampling Points for Spec. El. Method
np_u = 3;
np_v = 3;

