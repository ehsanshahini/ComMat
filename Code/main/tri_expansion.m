function aa=tri_expansion(dia)

toe=toeplitz(0:dia);
ind_u=triu(true(dia+1,dia+1));
ind_l=tril(true(dia+1,dia+1));
toe2=rot90(toe.*ind_l);
ind_l=rot90(ind_l);
aa=repmat(toe(ind_u),1,3)';
bb=repmat(toe2(ind_l),1,3)';
aa=[aa(:) bb(:)];
end
