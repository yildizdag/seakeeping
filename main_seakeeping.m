%==========================================================================
% Method        : Boundary Element Method (BEM) for Linear Seakeeping
% Discretization: Panel Method (Constant / Linear / Higher-Order Elements)
% Model         : Linear Potential Flow, Zero Forward Speed, Freq. Domain
% Green Function: Free-Surface Green Function (Telste & Noblesse)
%
% Description:
%  This code performs frequency-domain radiation and diffraction analyses
%  for floating objects using a direct BEM formulation. The solver computes
%  influence matrices based on the free-surface Green function and obtaines 
%  added mass, radiation damping, wave exciting forces, and rigid-body RAOs 
%  for the zero-speed seakeeping case.
%
% Author        : M. Erden Yildizdag
% Affiliation   : Istanbul Technical University
%========================================================================== 
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'cylinder_';
numPatch = 3;
%Parameters:
rho = 1000;
%-DOF per Node:
local_dof = 1;
%-Order:
p = 0; % 0:constant / 1:linear / 2:quadratic
%----------------
% Pre-Processing
%----------------
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
% iga2DmeshPlotNURBS(Nurbs2D);
%
bem2D = bem2Dmesh(Nurbs2D,p,local_dof);
bem2D.rho = rho;
bem2D.p = p; %-Element Order
bem2D.cg = [0.0,0.0,0.0]; %-Center of Gravity
bem2D.Lchar = 1;
bem2D.beta = 180;
bem2D.zeta = 0.5;
%----------
% Solution
%----------
%-Define the frequencies:
freq = [0.1,0.5:1:9.5,11:1.5:26,27.6409];
L = 1;
for i = 1:length(freq)
    disp(i)
    %-BEM Matrices:
    [H,G,C,bn] = bem_seakeeping(bem2D,freq(i),L);
    phi = (C-H)\(G*bn);
    bem2D.phi = phi;
    bem2D.bn = bn;
    [A,B] = addedMass_seakeeping(bem2D);
    fileAdded = fopen([FileName,num2str(freq(i)),'_added.txt'],'w');
    fileDamping = fopen([FileName,num2str(freq(i)),'_damping.txt'],'w');
 	fprintf(fileAdded,'%6.2f\n',freq(i));
	fprintf(fileDamping,'%6.2f\n',freq(i));
    for k = 1:6
        for l = 1:6
            fprintf(fileAdded,'%12.8f\t',A(k,l));
            fprintf(fileDamping,'%12.8f\t',B(k,l));
        end
        fprintf(fileAdded,'\n');
        fprintf(fileDamping,'\n');
    end
    fclose(fileAdded);
    fclose(fileDamping);
end
%-----------------
% Post-Processing
%-----------------
L = 1.284;
g = 9.81;
R = 0.18;
rho = 1000;
Adata = zeros(length(freq),6);
Bdata = zeros(length(freq),6);
for i = 1:length(freq)
    file1 = [FileName,num2str(freq(i)),'_added.txt'];
    file2 = [FileName,num2str(freq(i)),'_damping.txt'];
    fileID1 = fopen(file1);
    fileID2 = fopen(file2);
    Af = fscanf(fileID1,'%f');
    Ad = reshape(diag(reshape(Af(2:end),6,6)),1,6);
    Adata(i,:) = Ad./(L*(R^3/2)*rho*pi^2);
    Bf = fscanf(fileID2,'%f');
    Bd = reshape(diag(reshape(Bf(2:end),6,6)),1,6);
    Bdata(i,:) = Bd.*freq(i)./(L*(R^3/2)*rho*pi^2).*sqrt(R/g);
end
f = freq.*sqrt(L/g);
figure;
subplot(1,2,1)
plot(f,Adata(:,3));
xlim([0 10])
subplot(1,2,2)
plot(f,Bdata(:,3));
xlim([0 10])