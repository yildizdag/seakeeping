function bem2D = bem2Dmesh(Nurbs2D,N,local_dof)
nel = 0;
for k = 1:Nurbs2D.numpatch
    nel = nel + Nurbs2D.nel{k};
end
bem2D.nel = nel;
bem2D.N = N;
bem2D.local_dof = local_dof;
%
nodeData = zeros(N*N*nel,3);
JacMatData = zeros(2,2,N*N*nel);
InvJacMatData = zeros(2,2,N*N*nel);
JacobianData = zeros(N*N*nel,1);
count_el = 1;
count_node = 1;
%
if N == 0
    xi = 0; eta = 0;
else
    xi = linspace(-1,1,N+1);
    eta = linspace(-1,1,N+1);
end
%
epsilon = 1E-4;
%
for k = 1:Nurbs2D.numpatch
    for el = 1:Nurbs2D.nel{k}
        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);
        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);
        u_sample = (0.5.*(1-xi).*u1+0.5.*(1+xi).*u2);
        v_sample = (0.5.*(1-eta).*v1+0.5.*(1+eta).*v2);
        count = 1;
        CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
        du = (Nurbs2D.knots.U{k}(iu+1)-Nurbs2D.knots.U{k}(iu))/2;
        dv = (Nurbs2D.knots.V{k}(iv+1)-Nurbs2D.knots.V{k}(iv))/2;
        for i = 1:N+1
            for j = 1:N+1
                dNu = dersbasisfuns(iu,u_sample(i),Nurbs2D.order{k}(1)-1,1,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(j),Nurbs2D.order{k}(2)-1,1,Nurbs2D.knots.V{k});
                [~,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{k}(1),Nurbs2D.order{k}(2),CP,1,1);
                nodeData(count_node,:) = epsilon.*(dS(:,1,1)'./epsilon);
                %
                a = du * dS(1,2,1);
                b = du * dS(2,2,1);
                c = dv * dS(1,1,2);
                d = dv * dS(2,1,2);
                %
                detJ = a*d - b*c;
                %
                JacMatData(:,:,count_node) = [a b; c d];
                JacobianData(count_node)   = detJ;
                InvJacMatData(:,:,count_node) = (1/detJ).*[d -b; -c a];
                count = count+1;
                count_node = count_node+1;
            end
        end
        count_el = count_el+1;
    end
end
%
TOL = 1e-4;
[nodes_sem, IA, IC] = uniquetol(nodeData, TOL, 'ByRows', true);
Jmat = JacMatData(:,:,IA);
InvJmat = InvJacMatData(:,:,IA);
J = JacobianData(IA);
elemNode = reshape(IC, (N+1)*(N+1), nel).';
conn_sem = zeros(nel, local_dof*(N+1)*(N+1));
for d = 1:local_dof
    conn_sem(:, d:local_dof:end) = local_dof*elemNode - (local_dof - d);
end
bem2D.nodes = nodes_sem;
bem2D.conn = conn_sem;
bem2D.Jmat = Jmat;
bem2D.J = J;
bem2D.InvJmat = InvJmat;