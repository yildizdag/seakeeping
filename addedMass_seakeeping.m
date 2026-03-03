function [A,B] = addedMass_seakeeping(bem2D)
%
A = zeros(6,6);
B = zeros(6,6);
[xgp,wgp,ngp] = gaussQuad2d(2,2);
[~, dN] = shapefunc2D(xgp(:,1),xgp(:,2),bem2D.p);
for i = 1:6
    if bem2D.p == 0
        for el = 1:bem2D.nel
            %
            phi_el = (bem2D.phi(el,i));
            %
            bn_el = bem2D.bn(el,1:6);
            %
            for gp = 1:ngp
                a1j = dN(1,:,gp)*bem2D.nodes(bem2D.conn(el,:),:);
                a2j = dN(2,:,gp)*bem2D.nodes(bem2D.conn(el,:),:);
                J = norm(cross(a1j,a2j));
                A(i,:) = A(i,:) + (1000*wgp(gp)*J)*(real(phi_el)).*bn_el;
                B(i,:) = B(i,:) + (1000*wgp(gp)*J)*(imag(phi_el)).*bn_el;
            end
        end
    end
end