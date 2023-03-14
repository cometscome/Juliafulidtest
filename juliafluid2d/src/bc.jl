function reflect_bc(f::Fluid2d{Ni,Nj,Nb,Nf},bc_type) where {Ni,Nj,Nb,Nf}
    if bc_type == "periodical_in_i"
		bc_periodical_in_i(f);
	else 
		error("ERROR Fluid2d::reflect_bc > BC-type cannnot be specified.")
    end
end

function bc_periodical_in_i(f::Fluid2d{Ni,Nj,Nb,Nf}) where {Ni,Nj,Nb,Nf}
    u = get_u(f)
    v = get_v(f)
    e = get_e(f)
    rho = get_rho(f)

    for i=1:Nb#(int i = 0; i < Nb; ++i) {
		for j=1:Nj#(int j = 0; j < Nj; ++j) {
			rho[j,i] = rho[j,Ni-2*Nb+i];
			u[j,i] = u[j,Ni-2*Nb+i];
			v[j,i] = v[j,Ni-2*Nb+i];
			e[j,i] = e[j,Ni-2*Nb+i];
        end
	end
	for i=1:Nb#(int i = 0; i < Nb; ++i) {
		for j=1:Nj#(int j = 0; j < Nj; ++j) {
			rho[j,Ni-Nb+i] = rho[j,Nb+i];
			u[j,Ni-Nb+i] = u[j,Nb+i];
			v[j,Ni-Nb+i] = v[j,Nb+i];
			e[j,Ni-Nb+i] = e[j,Nb+i];
        end
	end
end