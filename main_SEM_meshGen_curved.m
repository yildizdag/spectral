%---------------------------------------------------------------
% NURBS-Enhanced Mesh Generator for
% Spectral Element Method (SEM)
%---------------------------------------------------------------
%
clc; clear; close all;
addpath('geometry')
% Read the Geometry imported from Rhino:
FileName = 'stiffner_';
numPatch = 5; %Enter # Patches
% Degrees of Freedom per each node:
local_dof = 1;
%-----------------------------------------------------------------
% Create 2-D Nurbs Structure (reads FileName)
%-----------------------------------------------------------------
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
%--------------------------------------
% Refinement (if necessary)
%--------------------------------------
ur = 0; % Refinement Level in u direction
vr = 0; % Refinement Level in v direction
Nurbs2D = hrefine2D(Nurbs2D,1,ur,vr);
Nurbs2D = hrefine2D(Nurbs2D,2,ur,vr);
Nurbs2D = hrefine2D(Nurbs2D,3,ur,vr);
Nurbs2D = hrefine2D(Nurbs2D,4,ur,vr);
Nurbs2D = hrefine2D(Nurbs2D,5,ur,vr);
%-----------------------------------------------------
% Plot Imported 2-D NURBS Structure
%-----------------------------------------------------
figure;
iga2DmeshPlotNURBS(Nurbs2D);
axis off
%------------------------------------------------------------------------------------
% Points for Spectral Element Method (e.g. 5 x 5, 3 x 3, etc.)
%------------------------------------------------------------------------------------
np_u = 5;
np_v = 5;
tot_el = 0;
for k = 1:Nurbs2D.numpatch
    tot_el = tot_el + Nurbs2D.nel{k};
end
elData = zeros(np_u*np_v,3,tot_el);
count_el = 0;
count_node = 1;
figure;
iga2DmeshPlotNURBS(Nurbs2D);
for k = 1:Nurbs2D.numpatch
    for el = 1:Nurbs2D.nel{k}
        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
        u_sample = linspace(Nurbs2D.knots.U{k}(iu),Nurbs2D.knots.U{k}(iu+1),np_u);
        v_sample = linspace(Nurbs2D.knots.V{k}(iv),Nurbs2D.knots.V{k}(iv+1),np_v);
        x_sample = zeros(np_u*np_v,1);
        y_sample = zeros(np_u*np_v,1);
        z_sample = zeros(np_u*np_v,1);
        count = 1;
        count_el = count_el+1;
        for j = 1:np_v
            for r = 1:np_u
                dNu = dersbasisfuns(iu,u_sample(r),Nurbs2D.order{k}(1)-1,0,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(j),Nurbs2D.order{k}(2)-1,0,Nurbs2D.knots.V{k});
                CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
                Sw = zeros(4,1);
                for i = 1:4
                    Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{k}(1),Nurbs2D.order{k}(2))*dNv(1,:)';
                end
                S = Sw(1:3,:) / Sw(4);
                x_sample(count) = S(1);
                y_sample(count) = S(2);
                z_sample(count) = S(3);
                elData(count,:,count_el) = [S(1) S(2) S(3)];
                count = count+1;
                count_node = count_node+1;
            end
        end
    end
end
%---------------------------------------------------------------
% Patch Connectivity
% Nodal Coordinates (nodes)
% Connectivity (conn)
%---------------------------------------------------------------
c_elData = zeros(4,3,count_el);
for i =1:count_el
    c_elData(:,:,i) = [elData(1,:,i);
                             elData(np_u,:,i);
                             elData(np_u*np_v,:,i);
                             elData(np_u*(np_v-1)+1,:,i)];
end
for i = 1:count_el-1
    el = 2;
    for j = el:count_el
        [~,ia,ib] = intersect(c_elData(:,:,i),c_elData(:,:,j),'rows');
        if length(ia) == 2
            if (ismember(1,ia) && ismember(2,ia))
                if (ismember(1,ib) && ismember(2,ib))
                    elData(1:np_u,:,j) = elData(1:np_u,:,i);
                elseif (ismember(2,ib) && ismember(3,ib))
                    elData(np_u:np_u:end,:,j) = elData(1:np_u,:,i);
                elseif (ismember(3,ib) && ismember(4,ib))
                    elData(np_u*(np_v-1)+1:np_u*np_v,:,j) = elData(1:np_u,:,i);
                elseif (ismember(1,ib) && ismember(4,ib))
                    elData(1:np_u:end,:,j) = elData(1:np_u,:,i);
                end
            elseif (ismember(2,ia) && ismember(3,ia))
                if (ismember(1,ib) && ismember(2,ib))
                    elData(1:np_u,:,j) = elData(np_u:np_u:end,:,i);
                elseif (ismember(2,ib) && ismember(3,ib))
                    elData(np_u:np_u:end,:,j) = elData(np_u:np_u:end,:,i);
                elseif (ismember(3,ib) && ismember(4,ib))
                    elData(np_u*(np_v-1)+1:np_u*np_v,:,j) = elData(np_u:np_u:end,:,i);
                elseif (ismember(1,ib) && ismember(4,ib))
                    elData(1:np_u:end,:,j) = elData(np_u:np_u:end,:,i);
                end
            elseif (ismember(3,ia) && ismember(4,ia))
                if (ismember(1,ib) && ismember(2,ib))
                    elData(1:np_u,:,j) = elData(np_u*(np_v-1)+1:np_u*np_v,:,i);
                elseif (ismember(2,ib) && ismember(3,ib))
                    elData(np_u:np_u:end,:,j) = elData(np_u*(np_v-1)+1:np_u*np_v,:,i);
                elseif (ismember(3,ib) && ismember(4,ib))
                    elData(np_u*(np_v-1)+1:np_u*np_v,:,j) = elData(np_u*(np_v-1)+1:np_u*np_v,:,i);
                elseif (ismember(1,ib) && ismember(4,ib))
                    elData(1:np_u:end,:,j) = elData(np_u*(np_v-1)+1:np_u*np_v,:,i);
                end
            elseif (ismember(1,ia) && ismember(4,ia))
                if (ismember(1,ib) && ismember(2,ib))
                    elData(1:np_u,:,j) = elData(1:np_u:end,:,i);
                elseif (ismember(2,ib) && ismember(3,ib))
                    elData(np_u:np_u:end,:,j) = elData(1:np_u:end,:,i);
                elseif (ismember(3,ib) && ismember(4,ib))
                    elData(np_u*(np_v-1)+1:np_u*np_v,:,j) = elData(1:np_u:end,:,i);
                elseif (ismember(1,ib) && ismember(4,ib))
                    elData(1:np_u:end,:,j) = elData(1:np_u:end,:,i);
                end
            end
        end
    end
    el = el + 1;
end
nodeData = reshape(permute(elData,[1,3,2]),[],size(elData,2));
TOL = 0.001;
nodes = uniquetol(nodeData,TOL,'ByRows',true);
conn = zeros(tot_el,np_u*np_v);
for i = 1:tot_el
    for j = 1:np_u*np_v
        node_id = find(ismembertol(nodes, elData(j,:,i),TOL,'ByRows',true));
        conn(i,j) = node_id;
    end
end
hold on
scatter3(nodes(:,1),nodes(:,2),nodes(:,3),80,'b','filled');
hold off