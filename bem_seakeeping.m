function [H,G,C,b] = bem_seakeeping(bem2D,we,L)
%
if bem2D.p == 0
    ns = bem2D.nel;
else
    ns = size(bem2D.nodes,1);
end
%-Initialize BEM Matrices:
H = zeros(ns,ns);
G = zeros(ns,ns);
C = zeros(ns,ns);
b = zeros(ns,ns);
%
if bem2D.p == 0
    for i = 1:ns
        xs = bem2D.snodes(i,:);
        for el = 1:bem2D.nel
            xf = bem2D.snodes(el,:);
            dist = norm(xf-xs)/sqrt(bem2D.A(i));
            if dist == 0
                [xgp,wgp,ngp] = gaussQuad2d(8,8);
            elseif dist < 1.5 && dist > 0
                [xgp,wgp,ngp] = gaussQuad2d(8,8);
            elseif dist < 3 && dist >= 1.5
                [xgp,wgp,ngp] = gaussQuad2d(6,6);
            elseif dist < 4.5 && dist >= 3
                [xgp,wgp,ngp] = gaussQuad2d(4,4);
            else
                [xgp,wgp,ngp] = gaussQuad2d(2,2);
            end
            %
            [N, dN] = shapefunc2D(xgp(:,1),xgp(:,2),bem2D.p);
            %
            if dist < 1E-3
                %-Sub elements in parametric domain:
                xi1 = [-1,-1;0,-1;-1,0;0,0];
                xi2 = [0,-1;1,-1;0,0;1,0];
                xi3 = [-1,0;0,0;-1,1;0,1];
                xi4 = [0,0;1,0;0,1;1,1];
                [N1,~] = shapefunc2D(xi1(:,1),xi1(:,2),bem2D.p);
                [N2,~] = shapefunc2D(xi2(:,1),xi2(:,2),bem2D.p);
                [N3,~] = shapefunc2D(xi3(:,1),xi3(:,2),bem2D.p);
                [N4,~] = shapefunc2D(xi4(:,1),xi4(:,2),bem2D.p);
                el_nodes1 = N1*bem2D.nodes(bem2D.conn(el,:),:);
                el_nodes2 = N2*bem2D.nodes(bem2D.conn(el,:),:);
                el_nodes3 = N3*bem2D.nodes(bem2D.conn(el,:),:);
                el_nodes4 = N4*bem2D.nodes(bem2D.conn(el,:),:);
                for gp = 1:ngp
                    xf1_gp = N(gp,:)*el_nodes1;
                    xf2_gp = N(gp,:)*el_nodes2;
                    xf3_gp = N(gp,:)*el_nodes3;
                    xf4_gp = N(gp,:)*el_nodes4;
                    a1j1 = dN(1,:,gp)*el_nodes1; a2j1 = dN(2,:,gp)*el_nodes1;
                    a1j2 = dN(1,:,gp)*el_nodes2; a2j2 = dN(2,:,gp)*el_nodes2;
                    a1j3 = dN(1,:,gp)*el_nodes3; a2j3 = dN(2,:,gp)*el_nodes3;
                    a1j4 = dN(1,:,gp)*el_nodes4; a2j4 = dN(2,:,gp)*el_nodes4;
                    J1 = norm(cross(a1j1,a2j1));
                    J2 = norm(cross(a1j2,a2j2));
                    J3 = norm(cross(a1j3,a2j3));
                    J4 = norm(cross(a1j4,a2j4));
                    %BEM Matrices:
                    Xdless = xs./(L); Xidless1 = xf1_gp./(L); Xidless2 = xf2_gp./(L); Xidless3 = xf3_gp./(L); Xidless4 = xf4_gp./(L);
                    f = we^2*L/9.81;
                    [GFS1,~] = freeSurfGreenFuncModern(Xidless1(1),Xidless1(2),Xidless1(3),Xdless(1),Xdless(2),Xdless(3),f);
                    [GFS2,~] = freeSurfGreenFuncModern(Xidless2(1),Xidless2(2),Xidless2(3),Xdless(1),Xdless(2),Xdless(3),f);
                    [GFS3,~] = freeSurfGreenFuncModern(Xidless3(1),Xidless3(2),Xidless3(3),Xdless(1),Xdless(2),Xdless(3),f);
                    [GFS4,~] = freeSurfGreenFuncModern(Xidless4(1),Xidless4(2),Xidless4(3),Xdless(1),Xdless(2),Xdless(3),f);
                    G(i,el) = G(i,el) + (-1*GFS1*wgp(gp)*J1) + (-1*GFS2*wgp(gp)*J2) + (-1*GFS3*wgp(gp)*J3) + (-1*GFS4*wgp(gp)*J4);
                end
            else
                for gp = 1:ngp
                    xf_gp = N(gp,:)*bem2D.nodes(bem2D.conn(el,:),:);
                    a1j = dN(1,:,gp)*bem2D.nodes(bem2D.conn(el,:),:);
                    a2j = dN(2,:,gp)*bem2D.nodes(bem2D.conn(el,:),:);
                    J = norm(cross(a1j,a2j));
                    nf = bem2D.n(el,:);
                    %BEM Matrices:
                    Xdless = xs./(L); Xidless = xf_gp./(L);
                    f = we^2*L/9.81;
                    [GFS,GxFS] = freeSurfGreenFuncModern(Xidless(1),Xidless(2),Xidless(3),Xdless(1),Xdless(2),Xdless(3),f);
                    H(i,el) = H(i,el) + (-1*(GxFS(1)*nf(1)+GxFS(2)*nf(2)+GxFS(3)*nf(3))*wgp(gp)*J);
                    G(i,el) = G(i,el) + (-1*GFS*wgp(gp)*J);
                end
            end
        end
    end
end
disp(det(G))