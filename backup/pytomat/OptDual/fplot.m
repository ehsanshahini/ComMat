function T=fplot(A,X,T,M)
% T=fplot(A,X,T,M)
% plot 3D patched trivalent molecules
% A is the connectivity matrix (sparse logical), 
% X the cartesian coordinates matrix (Natoms-by-3)
% M specifies the largest ring being searched, default is set to 8.

if nargin<3
    T=cell(6,1);
end
if nargin<4
    M=8;
end
if nargin>4
    error 'too many input argument!'
end

% figure
set(gcf,'Units','normalized')
% set(gcf,'position',[0 0 1 1])
% set(gcf,'position',[1 0 1 1])
set(gcf,'numbertitle','off','name',['C' num2str(length(X))]);
gplot3(A,X,'-','LineWidth',1,'color',[1 1 1]*0.1); axis equal;axis off
FA=0.8;

if isempty(T)||nargin<3
    T=ring_ind_c(A,M);
end
for k=[1:3 5:M-2]
    if numel(T{k}),patch('Faces',T{k},'Vertices',X,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',FA);end
end
if numel(T{4}),patch('Faces',T{4},'Vertices',X,'FaceColor',[1 1 1],'FaceAlpha',FA);end

axis vis3d
rotate3d on

if nargout<1
    clear T
end
end
