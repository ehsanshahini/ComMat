function tileF(n)
% tileF(n)
% this function rearranges all the currently existing figures so that they
% appear like rectangular grids of plots
%
% n can be a two-component vector or a scalar. In the former case, this
% function tiles the figures in an n(1)-by-n(2) rectangular pattern that
% fills the screen. For n being a scalar, the situation is the same as
% requiring n(2)=n(1) in the former case, i.e. an n-by-n lattice.

% keyboard
fig=allchild(0);
nfig=length(fig);
err=1e-2;
if ~nargin
%     n=[1 1]*ceil(sqrt(nfig)); % square tiling
    n=[1 nfig]; % linear tiling
else
    if numel(n)==1
        n=[n n];
    end
    while n(1)*n(2)<nfig
        n(2)=n(2)+1;
    end
end

fig=fig(end:-1:1);
width=1/n(2);
if isunix % mac
    height=0.80125/n(1);
    base=0.07875;
else
    height=1/n(1);
    base=0;
end
for k=1:nfig
    set(fig(k),'Units','normalized')
%     set(fig(k),'Position',[floor(k/n(1)-err)/n(2) mod(k-1,n(1))/n(1) 1/n(2) 1/n(1)])
    set(fig(k),'Position',[floor(k/n(1)-err)*width mod(k-1,n(1))*height+base width height])
end
end