function [H,G,C,b] = bem_seakeeping(bem2D,we,L)
%
if bem2D.N == 0
    ns = bem2D.nel;
end
%-Initialize BEM Matrices:
H = zeros(ns,ns);
G = zeros(ns,ns);
C = zeros(ns,ns);
b = zeros(ns,ns);
%
count = 0;
%
if bem2D.N == 0
    for i = 1:ns
        count = count + 1;
        xs = bem2D.snodes(i,:);
        for el = 1:bem2D.nel
            xf = bem2D.snodes(el,:);
            dist = norm(xf-xs)/sqrt(bem2D.A(i));
            if dist == 0
                [xgp,wgp,ngp] = gaussQuad2d(6,6);
            elseif dist < 1 && dist > 0
                [xgp,wgp,ngp] = gaussQuad2d(6,6);
            elseif dist < 2 && dist >= 1
                [xgp,wgp,ngp] = gaussQuad2d(4,4);
            elseif dist < 0.15 && dist >= 0.1
                [xgp,wgp,ngp] = gaussQuad2d(4,4);
            else
                [xgp,wgp,ngp] = gaussQuad2d(2,2);
            end
            %
            [N, dN] = shapefunc2D(xgp(:,1),xgp(:,2),bem2D.N);
            %
            for gp = 1:ngp
                xf_gp = N(gp,:)*sem2D.nodes(sem2D.conn(el,:),:);
                a1j = dN(1,:,gp)*sem2D.nodes(sem2D.conn(el,:),:);
                a2j = dN(2,:,gp)*sem2D.nodes(sem2D.conn(el,:),:);
                J = norm(cross(a1j,a2j));
                %BEM Matrices:
                Xdless = xs./(L); Xidless = xf_gp./(L);
                f = we^2*L/9.81;
                [GFS,GxFS] = freeSurfGreenFuncModern(Xidless(1),Xidless(2),Xidless(3),Xdless(1),Xdless(2),Xdless(3),f);
                H(i,el) = H(i,el) + (-1*(GxFS(1)*n2(1)+GxFS(2)*n2(2)+GxFS(3)*n2(3))*wgp(gp)*J).*N(gp,:);
                G(i,el) = G(i,el) + (-1*GFS*wgp(gp)*J).*N(gp,:);
            end
        end
    end
end