function [A,B] = addedMass_cylindricalShell(Nurbs2D)
modeNum = Nurbs2D.modeNum;
rigidModeNum = Nurbs2D.rigidModeNum;
A = zeros(rigidModeNum+2*modeNum,rigidModeNum+2*modeNum);
B = zeros(rigidModeNum+2*modeNum,rigidModeNum+2*modeNum);
[xgp,wgp,ngp] = gaussQuad2d(3,3);
for i = 1:rigidModeNum+2*modeNum
    for kk = 3:5
        for el = 1:Nurbs2D.nel{kk}
            %Parametric Domain of the Element:
            iu = Nurbs2D.INC{kk}(Nurbs2D.IEN{kk}(1,el),1);   
            iv = Nurbs2D.INC{kk}(Nurbs2D.IEN{kk}(1,el),2);
            %Connectivity:
            lm_loc = reshape(fliplr(Nurbs2D.LM{kk}(:,:,el)),1,Nurbs2D.nen{kk});
            %Potential on the Element:
            phi_local = (Nurbs2D.phi(lm_loc,i));
            %Modal Displacement:
            eigvec_local = Nurbs2D.b(lm_loc,:);
            %Control Points:
            CP = Nurbs2D.cPoints{kk}(:,iu-Nurbs2D.order{kk}(1)+1:iu,iv-Nurbs2D.order{kk}(2)+1:iv);
            %Gaussian Quadrature over Element:
            for jj = 1:ngp
                %Parametric Location:
                cu1 = (Nurbs2D.knots.U{kk}(iu) + Nurbs2D.knots.U{kk}(iu+1));
                cu2 = (Nurbs2D.knots.U{kk}(iu+1) - Nurbs2D.knots.U{kk}(iu));
                u = (cu1+cu2*xgp(jj,1))/2;
                cv1 = (Nurbs2D.knots.V{kk}(iv) + Nurbs2D.knots.V{kk}(iv+1));
                cv2 = (Nurbs2D.knots.V{kk}(iv+1) - Nurbs2D.knots.V{kk}(iv));
                v = (cv1+cv2*xgp(jj,2))/2;
                %Shape Functions:
                dNu = dersbasisfuns(iu,u,Nurbs2D.order{kk}(1)-1,2,Nurbs2D.knots.U{kk});
                dNv = dersbasisfuns(iv,v,Nurbs2D.order{kk}(2)-1,2,Nurbs2D.knots.V{kk});
                [dR,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{kk}(1),Nurbs2D.order{kk}(2),CP,1,1);
                %From Physical to Parametric:
                a1 = dS(:,2,1);
                a2 = dS(:,1,2);
                J1 = norm(cross(a1,a2));
                %From parametric to parent
                J2 = diag([Nurbs2D.knots.U{kk}(iu+1) - Nurbs2D.knots.U{kk}(iu),...
                                     Nurbs2D.knots.V{kk}(iv+1) - Nurbs2D.knots.V{kk}(iv)])/2;
                % Jacobian
                J = J1*det(J2);
                %Reshape Shape Functions:
                dRB = reshape(dR(:,:,1,1),1,Nurbs2D.nen{kk});
                phi_local1 = dRB*real(phi_local);
                phi_local2 = dRB*imag(phi_local);
                eigvec_local1 = dRB*eigvec_local;
                %Added Mass Matrix:
                for j = 1:rigidModeNum+2*modeNum
                    A(j,i) = A(j,i) + (1000*wgp(jj)*J).*(phi_local1*(eigvec_local1(j)));
                    B(j,i) = B(j,i) + (1000*wgp(jj)*J).*(phi_local2*(eigvec_local1(j)));
                end
            end
        end
    end
end