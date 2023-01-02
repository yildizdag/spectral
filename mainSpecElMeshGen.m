%---------------------------------------------------------------
% NURBS-Enhanced Mesh Generator for
% Spectral Element Method
%---------------------------------------------------------------
%
% Read the Geometry imported from Rhino:
FileName = 'elliptical_cutout_';
numPatch = 6; %Enter # Patches
% Degrees of Freedom per each node:
local_dof = 1;
% CREATE 2D IGA MESH (reads FileName):
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
%---------------------------------------------------------------
% Plot Imported 2-D NURBS Structure
iga2DmeshPlotNURBS(Nurbs2D);
%---------------------------------------------------------------
% Points for Spec. El. Method
% 5 x 5, 3 x 3, etc.
np_u = 3;
np_v = 3;
tot_el = 0;
for k = 1:Nurbs2D.numpatch
    tot_el = tot_el + Nurbs2D.nel{k};
end
elData = zeros(np_u*np_v,3,tot_el);
nodeData = zeros(np_u*np_v*tot_el,3);
count_el = 1;
count_node = 1;
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
                nodeData(count_node,:) = [S(1), S(2), S(3)];
                count = count+1;
                count_node = count_node+1;
            end
        end
        hold on
        scatter3(x_sample,y_sample,z_sample,80,'b','filled');
        hold off
        count_el = count_el+1;
    end
end
%---------------------------------------------------------------
% Nodal Coordinates (nodes)
% Connectivity (conn)
%---------------------------------------------------------------
nodes = uniquetol(nodeData,'ByRows',true);
conn = zeros(tot_el,np_u*np_v);
for i = 1:tot_el
    for j = 1:np_u*np_v
        node_id = find(ismembertol(nodes, elData(j,:,i),1E-6,'ByRows',true));
        conn(i,j) = node_id;
    end
end