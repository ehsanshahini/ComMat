function [E,E_gradient,screw_gradient]=helix_Bead_E(R,ii,k,n)
% energy potential form proposed by BYJ lab, NTUCH et al. 2006 Nov. 24 18:56
% (t0,u0) are the screw symmetry parameters

t0=R(end-1);
u0=R(end);
R=reshape(R(1:end-2),[],3);
N=size(R,1);

e0 = 50;    %% parameters are account in eV
e1 = 9.3;
b0 = 1.43;

R=[
   bsxfun(@plus,R,[0 0 -u0])*RotatM([t0 0 0 1]);
   R;
   bsxfun(@plus,R,[0 0 u0])*RotatM([t0 0 0 -1]);
   ];

for kk=2:(n-1)/2
    R=[R;
       bsxfun(@plus,R(N+1:N*2,:),[0 0 -u0*kk])*RotatM([t0*kk 0 0 1]);
       bsxfun(@plus,R(N+1:N*2,:),[0 0 u0*kk])*RotatM([t0*kk 0 0 -1]);
      ];
end
% A=false(N*n,N*n);
% A(ii(:,2)+(ii(:,3)-1)*N*n)=true;
% A=reshape(A,N*n,N*n);
% A=A|A';
% fplot(Ak,R);
% return

uu=R(ii(:,2),:)-R(ii(:,3),:);
u1=uu(k(1:N),:);
u2=uu(k(N+1:2*N),:);
u3=uu(k(2*N+1:3*N),:);

B1 = repmat(sqrt(sum(u1.^2,2)),1,3); % bond length
B2 = repmat(sqrt(sum(u2.^2,2)),1,3);
B3 = repmat(sqrt(sum(u3.^2,2)),1,3);

uu=bsxfun(@rdivide,uu,sqrt(sum(uu.^2,2)));
    
u1 = u1./B1; % (bonding) unit vectors
u2 = u2./B2;
u3 = u3./B3;

stretching = 0.5*e0*sum((B1(:,1)-b0).^2+(B2(:,1)-b0).^2+(B3(:,1)-b0).^2);
bending = (u1+u2).^2+(u2+u3).^2+(u3+u1).^2; bending=0.5*e1*sum(bending(:));

E = stretching + bending;

if nargout > 1
    %% bond stretching term
    E0_gradient=u1.*(B1-b0)+u2.*(B2-b0)+u3.*(B3-b0);

    u11=uu(k(3*N+1:4*N),:);u12=uu(k(4*N+1:5*N),:);
    u21=uu(k(5*N+1:6*N),:);u22=uu(k(6*N+1:7*N),:);
    u31=uu(k(7*N+1:8*N),:);u32=uu(k(8*N+1:9*N),:);
%     u11=u11./repmat(sqrt(sum(u11.^2,2)),1,3);u12=u12./repmat(sqrt(sum(u12.^2,2)),1,3);
%     u21=u21./repmat(sqrt(sum(u21.^2,2)),1,3);u22=u22./repmat(sqrt(sum(u22.^2,2)),1,3);
%     u31=u31./repmat(sqrt(sum(u31.^2,2)),1,3);u32=u32./repmat(sqrt(sum(u32.^2,2)),1,3);
    
    u1pu2=u1+u2;u2pu3=u2+u3;u3pu1=u3+u1;
    u1xty=u1.*u1(:,[2 3 1])./B1;
    u1ztx=u1(:,[3 1 2]).*u1./B1;
    u2xty=u2.*u2(:,[2 3 1])./B2;
    u2ztx=u2(:,[3 1 2]).*u2./B2;
    u3xty=u3.*u3(:,[2 3 1])./B3;
    u3ztx=u3(:,[3 1 2]).*u3./B3;
    u1res=(1-u1.^2)./B1;
    u2res=(1-u2.^2)./B2;
    u3res=(1-u3.^2)./B3;
    u1diff=2*u1-u11-u12;
    u2diff=2*u2-u21-u22;
    u3diff=2*u3-u31-u32;
    %% bond bending term
    E1_gradient      = u1pu2.*(u1res+u2res) ...
                     - u1pu2(:,[2 3 1]).*(u1xty+u2xty) ...
                     - u1pu2(:,[3 1 2]).*(u1ztx+u2ztx) ...
                     + u2pu3.*(u2res+u3res) ...
                     - u2pu3(:,[2 3 1]).*(u2xty+u3xty) ...
                     - u2pu3(:,[3 1 2]).*(u2ztx+u3ztx) ...
                     + u3pu1.*(u3res+u1res) ...
                     - u3pu1(:,[2 3 1]).*(u3xty+u1xty) ...
                     - u3pu1(:,[3 1 2]).*(u3ztx+u1ztx) ...
                     + u1diff.*u1res ...
                     - u1diff(:,[2 3 1]).*u1xty ...
                     - u1diff(:,[3 1 2]).*u1ztx ...
                     + u2diff.*u2res ...
                     - u2diff(:,[2 3 1]).*u2xty ...
                     - u2diff(:,[3 1 2]).*u2ztx ...
                     + u3diff.*u3res ...
                     - u3diff(:,[2 3 1]).*u3xty ...
                     - u3diff(:,[3 1 2]).*u3ztx;

                        
    E_gradient = 2*e0*E0_gradient + e1*E1_gradient;
    
    E_gradient = -E_gradient;
    screw_gradient=[0;0];
    dx=1e-8;
    screw_gradient(1)=(helix_Bead_E([reshape(R(1:N,:),[],1);t0*(dx+1);u0],ii,k,n)-E)/dx/t0;
    screw_gradient(2)=(helix_Bead_E([reshape(R(1:N,:),[],1);t0;u0*(dx+1)],ii,k,n)-E)/dx/u0;
end
end
