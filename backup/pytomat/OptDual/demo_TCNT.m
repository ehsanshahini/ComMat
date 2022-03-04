
% this script randomly chooses a TCNT to generate and performs optimization
close all
clear
Nmax=100; % maximum number of atoms per rotational unit cell

% randomizing the carbon nanotorus considered
type=randperm(2);type=type(1); % Dnd or Dnh symmetry
% ncell=randperm(5)+3;ncell=ncell(1); % rotational symmetry number
ncell=5; % rotational symmetry number
switch type
    case 1
        ind=generate_Dnd_list(ncell,Nmax,0);
        ind=ind(randperm(size(ind,1)),:);ind=ind(1,:);
        sym='d';
        [R,A]=Dnd_TCNT(ind(2:end),ncell);
    case 2
        ind=generate_Dnh_list(ncell,Nmax,0);
        ind=ind(randperm(size(ind,1)),:);ind=ind(1,:);
        sym='h';
        [R,A]=Dnh_TCNT(ind(2:end),ncell);
end
M=8;
T=ring_ind_c(A,M); % ring indices of the real space structure

[Rd,Ad,Nd,Td]=Real2Dual(R,T); % get the dual space structure

tic
Xd = BFGS_Dual(Rd,Ad,1,@(R,n)R); % geometry optimization in dual space
Rn_Duall = squeeze(mean(reshape(Xd(Td',:),[3 size(R,1) 3]),1)); % map back to real space structure
t_d=toc;

tic
% Rn_VSEPR = BFGS_VSEPR(R,A,ncell,@sym_Cn); % making use of the rotational symmetry
Rn_VSEPR = BFGS_VSEPR(R,A,1,@sym_Cn); % not making use of the symmetry
t_r=toc;

figure
fplot(A,R,T);
title('Pre-optimized')
figure
fplot(A,Rn_Duall,T);
title(sprintf('Dual space optimized, took %f seconds',t_d))
figure
fplot(A,Rn_VSEPR,T);
title(sprintf('Real space optimized, took %f seconds',t_r))
tileF
