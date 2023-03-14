function input_coordinate!(f::Fluid2d{Ni,Nj,Nb,Nf},f_name) where {Ni,Nj,Nb,Nf}
    data = readlines(f_name)
    it = 0
    for i=1:Ni
        for j=1:Nj
            it += 1
            u = split(data[it])
            f.x[j,i] = parse(Float64,u[1])
        end
    end
    for i=1:Ni
        for j=1:Nj
            it += 1
            u = split(data[it])
            f.y[j,i] = parse(Float64,u[1])
        end
    end
    println("Fluid2d > Coordinate successfully read.")

end

function input_basic!(f::Fluid2d{Ni,Nj,Nb,Nf},f_name) where  {Ni,Nj,Nb,Nf}
    data = readlines(f_name)
    it = 0
    for i=1:Ni
        for j=1:Nj
            it += 1
            u = split(data[it])
            f.rho[j,i] = parse(Float64,u[1])
        end
    end
    for i=1:Ni
        for j=1:Nj
            it += 1
            u = split(data[it])
            f.u[j,i] = parse(Float64,u[1])
        end
    end
    for i=1:Ni
        for j=1:Nj
            it += 1
            u = split(data[it])
            f.v[j,i] = parse(Float64,u[1])
        end
    end
    for i=1:Ni
        for j=1:Nj
            it += 1
            u = split(data[it])
            f.e[j,i] = parse(Float64,u[1])
        end
    end
    println("Fluid2d > Basic variables successfully read.")
end

function calc_metrices_dx!(f::Fluid2d{Ni,Nj,Nb,Nf})  where  {Ni,Nj,Nb,Nf}
    for i=1:(Ni-2*Nb+2)#(int i = 0; i < Ni-2*Nb+2; ++i) {
		for j=1:(Nj-2*Nb+2)#(int j = 0; j < Nj-2*Nb+2; ++j) {
			f.ixS[j,i] = 0.5 * (f.y[Nb+j,Nb-1+i] - f.y[Nb-2+j,Nb-1+i])
			f.iyS[j,i] = -0.5 * (f.x[Nb+j,Nb-1+i] - f.x[Nb-2+j,Nb-1+i])
			f.jxS[j,i] = -0.5 * (f.y[Nb-1+j,Nb+i] - f.y[Nb-1+j,Nb-2+i])
			f.jyS[j,i] = 0.5 * (f.x[Nb-1+j,Nb+i] - f.x[Nb-1+j,Nb-2+i])
        end
	end

    for i=1:Ni-1#(int i = 0; i < Ni-1; ++i) {
		for j=1:Nj-1 #(int j = 0; j < Nj-1; ++j) {
			f.S[j,i] = ((f.x[1+j,1+i] - f.x[j,i]) * (f.y[1+j,i] - f.y[j,1+i]) -
							 (f.y[1+j,1+i] - f.y[j,i]) * (f.x[1+j,i] -f.x[j,1+i])) * 
							0.5
        end
	end

    for i=1:( Ni-2*Nb) #(int i = 0; i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb)#(int j = 0; j < Nj-2*Nb; ++j) {
			dx1 = f.x[Nb+j,Nb+1+i] - f.x[Nb+j,Nb+i]
			dy1 = f.y[Nb+j,Nb+1+i] - f.y[Nb+j,Nb+i]
			f.dx[j,i] = dx1 * dx1 + dy1 * dy1
			dx1 = f.x[Nb+1+j,Nb+i] - f.x[Nb+j,Nb+i]
			dy1 = f.y[Nb+1+j,Nb+i] - f.y[Nb+j,Nb+i]
			dd = dx1 * dx1 + dy1 * dy1
			f.dx[j,i] = min(f.dx[j,i], dd)
			dx1 = f.x[Nb+j,Nb-1+i] - f.x[Nb+j,Nb+i]
			dy1 = f.y[Nb+j,Nb-1+i] - f.y[Nb+j,Nb+i]
			dd = dx1 * dx1 + dy1 * dy1
			f.dx[j,i] = min(f.dx[j,i], dd)
			dx1 = f.x[Nb-1+j,Nb+i] - f.x[Nb+j,Nb+i];
			dy1 = f.y[Nb-1+j,Nb+i]- f.y[Nb+j,Nb+i];
			dd = dx1 * dx1 + dy1 * dy1;
			f.dx[j,i] = min(f.dx[j,i], dd);
			f.dx[j,i] = sqrt(f.dx[j,i]);
        end
	end


end

function output_basic(f::Fluid2d{Ni,Nj,Nb,Nf},tstep, 
    iter
    )  where  {Ni,Nj,Nb,Nf}
    filename = "b"*lpad(tstep,7,"0")*".dat"
    dir_o = get_dir_o(f)
    fp = open(dir_o * filename,"w")
    u = get_u(f)
    v = get_v(f)
    e = get_e(f)
    rho = get_rho(f)

    for i=1:Ni#(int i = 0;i < Ni;++i) {
		for j=1:Nj#(int j = 0;j < Nj;++j) {
			println(fp,rho[j,i]) # << "\n";
        end
	end
    for i=1:Ni#(int i = 0;i < Ni;++i) {
		for j=1:Nj#(int j = 0;j < Nj;++j) {
			println(fp,u[j,i]) # << "\n";
        end
	end
    for i=1:Ni#(int i = 0;i < Ni;++i) {
		for j=1:Nj#(int j = 0;j < Nj;++j) {
			println(fp,v[j,i]) # << "\n";
        end
	end
    for i=1:Ni#(int i = 0;i < Ni;++i) {
		for j=1:Nj#(int j = 0;j < Nj;++j) {
			println(fp,e[j,i]) # << "\n";
        end
	end
    println("tstep = $tstep, iter = $iter")

end