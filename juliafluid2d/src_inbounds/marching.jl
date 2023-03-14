

using .Fnd

function calc_cfl(f::Fluid2d{Ni,Nj,Nb,Nf},FL_coeff)  where  {Ni,Nj,Nb,Nf}
    rho = get_rho(f)
    u = get_u(f)
    v = get_v(f)
    e = get_e(f)
    dx = get_dx(f)


    nu = 1.0e10
	@inbounds for i=1:(Ni - 2 * Nb)#(int i = 0; i < Ni - 2 * Nb; ++i) {
		for j=1:( Nj - 2 * Nb)#(int j = 0; j < Nj - 2 * Nb; ++j) {
			tmp = dx[j,i] /
				 (calc_cs(f,rho[j + Nb,i + Nb], u[j + Nb,i + Nb],
					v[j + Nb,i + Nb], e[j + Nb,i + Nb]) +
					 sqrt(u[j + Nb,i + Nb] * u[j + Nb,i + Nb] + 
						+ v[j + Nb,i + Nb] * v[j + Nb,i + Nb]))
			nu = min(nu, tmp)
        end
	end
    return FL_coeff * nu
end

using InteractiveUtils

#=
function update_SA!(f::Fluid2d{Ni,Nj,Nb,Nf})where {Ni,Nj,Nb,Nf}
    S = get_S(f)
    u = get_u(f)
    v = get_v(f)
    e = get_e(f)
    rho = get_rho(f)
    arr_Q0= get_arr_Q0(f)

    for i=1:( Ni-2*Nb)#(int i = 0; i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb)#(int j = 0; j < Nj-2*Nb; ++j) {
            SA = 0.25 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j - 1,Nb + i] +
                S[Nb + j,Nb + i - 1] + S[Nb + j,Nb + i])
            Q0 = view(arr_Q0,:,j,i)
            calc_conserved!(rho[Nb + j,Nb + i], u[Nb + j,Nb + i], 
                v[Nb + j,Nb + i], e[Nb + j,Nb + i], SA, Q0)#view(arr_Q0,:,j,i))

        end
    end
end
=#

function march_ssprk3!(f::Fluid2d{Ni,Nj,Nb,Nf},dt,
        bc_type,reconstruction,flux_scheme) where {Ni,Nj,Nb,Nf}
    S = get_S(f)
    u = get_u(f)
    v = get_v(f)
    e = get_e(f)
    rho = get_rho(f)
    arr_Q0= get_arr_Q0(f)
    arr_Q1= get_arr_Q1(f)
    arr_Q2= get_arr_Q2(f)

    
    #update_SA!(f)

    
    @inbounds for i=1:( Ni-2*Nb)#(int i = 0; i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb)#(int j = 0; j < Nj-2*Nb; ++j) {
            SA = 0.25 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j - 1,Nb + i] +
                S[Nb + j,Nb + i - 1] + S[Nb + j,Nb + i])

            calc_conserved!(rho[Nb + j,Nb + i], u[Nb + j,Nb + i], 
            v[Nb + j,Nb + i], e[Nb + j,Nb + i], SA, view(arr_Q0,:,j,i))

        end
    end
    



    rhs(f,rho, u, v, e, arr_Q1, reconstruction, flux_scheme)
    
    @inbounds for i=1:(Ni-2*Nb)#(int i = 0;i < Ni-2*Nb; ++i) {
		for j=1:( Nj-2*Nb) #(int j = 0; j < Nj-2*Nb; ++j) {


			for k=1:Nf #(int k = 0; k < Nf; ++k) {
				arr_Q1[k,j,i] = arr_Q0[k,j,i] + dt * arr_Q1[k,j,i];
            end
			SA = 0.25 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j - 1,Nb + i] + 
				 S[Nb + j,Nb + i - 1] + S[Nb+ j,Nb + i])
            #println("$i $j $SA")

            rho[Nb + j,Nb + i],u[Nb + j,Nb + i], v[Nb+j,Nb + i], e[Nb+j,Nb + i] = calc_basic(f,view(arr_Q1,:,j,i), SA)
            #println((i,j,rho[Nb + j,Nb + i],u[Nb + j,Nb + i], v[Nb+j,Nb + i], e[Nb+j,Nb + i]))
            
		end
	end
    reflect_bc(f,bc_type)

    rhs(f,rho, u, v, e, arr_Q2, reconstruction, flux_scheme)

    @inbounds for i=1:(Ni-2*Nb) #(int i = 0;i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb) #(int j = 0; j < Nj-2*Nb; ++j) {
			for k=1:Nf#(int k = 0; k < Nf; ++k) {
				arr_Q2[k,j,i] = 0.75 * arr_Q0[k,j,i] + 
				  0.25 * (arr_Q1[k,j,i] + dt * arr_Q2[k,j,i]);
            end
            
			SA = 0.25 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j - 1,Nb + i] + 
				S[Nb + j,Nb + i - 1] + S[Nb + j,Nb + i])

            rho[Nb + j,Nb + i], 
				u[Nb + j,Nb + i], v[Nb + j,Nb + i], e[Nb+j,Nb + i]= calc_basic(f,view(arr_Q2,:,j,i) , SA, )

		end

	end
    reflect_bc(f,bc_type)

    rhs(f,rho, u, v, e, arr_Q1, reconstruction, flux_scheme);

    @inbounds for i=1:(Ni-2*Nb)#(int i = 0;i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb) #(int j = 0; j < Nj-2*Nb; ++j) {
			for k=1:Nf #(int k = 0; k < Nf; ++k) {
				arr_Q0[k,j,i] = arr_Q0[k,j,i] / 3.0 + 
				 2.0 / 3.0 * (arr_Q2[k,j,i]+ dt * arr_Q1[k,j,i])
            end
			SA = 0.25 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j - 1,Nb + i] + 
				S[Nb + j,Nb + i - 1] + S[Nb + j,Nb + i])
                rho[Nb + j,Nb + i], 
				u[Nb + j,Nb + i], v[Nb + j,Nb + i], e[Nb + j,Nb + i] = calc_basic(f,view(arr_Q0,:,j,i), SA, )
		end
	end
    
    reflect_bc(f,bc_type)
    add_t!(f,dt)
    #f.t += dt
    add_tstep!(f,1)
    #f.tstep += 1

end

function rhs(f::Fluid2d{Ni,Nj,Nb,Nf},rho,u,v,e,arr_Q,reconstruction,flux_scheme) where {Ni,Nj,Nb,Nf}
    S = get_S(f)
    ixS = get_ixS(f)
    iyS = get_iyS(f)
    jxS = get_jxS(f)
    jyS = get_jyS(f)
    #u = get_u(f)
    #v = get_v(f)
    #e = get_e(f)
    #rho = get_rho(f)
    arr_Fi = get_arr_Fi(f)
    arr_Fj = get_arr_Fj(f)

    QL = zeros(Nf)
    QR = zeros(Nf)
    FL = zeros(Nf)
    FR = zeros(Nf)
    LAM = zeros(Nf)
    R = zeros(Nf,Nf)
    Rinv = zeros(Nf,Nf)


    @inbounds for i=1:(Ni-2*Nb+1)#(int i = 0; i < Ni-2*Nb+1; ++i) {
		for j=1:(Nj-2*Nb)#(int j = 0; j < Nj-2*Nb; ++j) {
            SA = 0.5 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j,Nb + i - 1]);
            ixSA = 0.5 * (ixS[j + 1,i] + ixS[j + 1,i + 1]);
			iySA = 0.5 * (iyS[j + 1,i] + iyS[j + 1,i + 1]);
			jxSA = 0.5 * (jxS[j + 1,i] + jxS[j + 1,i + 1]);
			jySA = 0.5 * (jyS[j + 1,i] + jyS[j + 1,i + 1]);


			if reconstruction == "MUSCL_minmod_basic"
				rhoL, uL, vL, eL, rhoR, uR, vR, eR = reconstruction_basic_muscl3( 
					rho[Nb+j,Nb+i-2], rho[Nb+j,Nb+i-1], 
					rho[Nb+j,Nb+i], rho[Nb+j,Nb+i+1],
					u[Nb+j,Nb+i-2], u[Nb+j,Nb+i-1], 
					u[Nb+j,Nb+i], u[Nb+j,Nb+i+1],
					v[Nb+j,Nb+i-2], v[Nb+j,Nb+i-1], 
					v[Nb+j,Nb+i], v[Nb+j,Nb+i+1],
					e[Nb+j,Nb+i-2], e[Nb+j,Nb+i-1], 
					e[Nb+j,Nb+i], e[Nb+j,Nb+i+1])
			elseif reconstruction == "MP5_basic" 
				rhoL, uL, vL, eL, rhoR, uR, vR, eR = reconstruction_basic_mp5(
					rho[Nb+j,Nb+i-3], rho[Nb+j,Nb+i-2],
					rho[Nb+j,Nb+i-1], rho[Nb+j,Nb+i],
					rho[Nb+j,Nb+i+1], rho[Nb+j,Nb+i+2],
					u[Nb+j,Nb+i-3], u[Nb+j,Nb+i-2],
					u[Nb+j,Nb+i-1], u[Nb+j,Nb+i],
					u[Nb+j,Nb+i+1], u[Nb+j,Nb+i+2],
					v[Nb+j,Nb+i-3], v[Nb+j,Nb+i-2],
					v[Nb+j,Nb+i-1], v[Nb+j,Nb+i],
					v[Nb+j,Nb+i+1], v[Nb+j,Nb+i+2],
					e[Nb+j,Nb+i-3], e[Nb+j,Nb+i-2],
					e[Nb+j,Nb+i-1], e[Nb+j,Nb+i],
					e[Nb+j,Nb+i+1], e[Nb+j,Nb+i+2])
			
			else 
                error("ERROR Fluid2d::rhs > Reconstruction scheme cannnot be specified.")
            end



            if flux_scheme == "Roe_FDS"
				roe_fds!(f,rhoL, uL, vL, eL, rhoR, uR, vR, eR,
					ixSA, iySA, SA, view(arr_Fi,:,j,i),
                    QL,QR,FL,FR,LAM,R,Rinv);
			
			else 
				error("ERROR Fluid2d::rhs > Flux scheme cannnot be specified.")
            end
        end
    end
    


	@inbounds for i=1:(Ni-2*Nb)#(int i = 0; i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb+1)#(int j = 0; j < Nj-2*Nb+1; ++j) {

			SA = 0.5 * (S[Nb + j - 1,Nb + i - 1] + S[Nb + j - 1,Nb + i]);

			ixSA = 0.5 * (ixS[j,i + 1] + ixS[j+1,i + 1]);
			iySA = 0.5 * (iyS[j,i + 1] + iyS[j+1,i + 1]);
			jxSA = 0.5 * (jxS[j,i + 1] + jxS[j+1,i + 1]);
			jySA = 0.5 * (jyS[j,i + 1] + jyS[j+1,i + 1]);

		
			if reconstruction == "MUSCL_minmod_basic" 
				rhoL, uL, vL, eL, rhoR, uR, vR, eR = reconstruction_basic_muscl3( 
					rho[Nb+j-2,Nb+i], rho[Nb+j-1,Nb+i], 
					rho[Nb+j,Nb+i], rho[Nb+j+1,Nb+i],
					u[Nb+j-2,Nb+i], u[Nb+j-1,Nb+i], 
					u[Nb+j,Nb+i], u[Nb+j+1,Nb+i],
					v[Nb+j-2,Nb+i], v[Nb+j-1,Nb+i], 
					v[Nb+j,Nb+i], v[Nb+j+1,Nb+i],
					e[Nb+j-2,Nb+i], e[Nb+j-1,Nb+i], 
					e[Nb+j,Nb+i], e[Nb+j+1,Nb+i]);
			
			elseif reconstruction == "MP5_basic"
				rhoL, uL, vL, eL, rhoR, uR, vR, eR = reconstruction_basic_mp5(
					rho[Nb+j-3,Nb+i], rho[Nb+j-2,Nb+i],
					rho[Nb+j-1,Nb+i], rho[Nb+j,Nb+i],
					rho[Nb+j+1,Nb+i], rho[Nb+j+2,Nb+i],
					u[Nb+j-3,Nb+i], u[Nb+j-2,Nb+i],
					u[Nb+j-1,Nb+i], u[Nb+j,Nb+i],
					u[Nb+j+1,Nb+i], u[Nb+j+2,Nb+i],
					v[Nb+j-3,Nb+i], v[Nb+j-2,Nb+i],
					v[Nb+j-1,Nb+i], v[Nb+j,Nb+i],
					v[Nb+j+1,Nb+i], v[Nb+j+2,Nb+i],
					e[Nb+j-3,Nb+i], e[Nb+j-2,Nb+i],
					e[Nb+j-1,Nb+i], e[Nb+j,Nb+i],
					e[Nb+j+1,Nb+i], e[Nb+j+2,Nb+i])
			
			else 
				error("ERROR Fluid2d::rhs > Reconstruction scheme cannnot be specified.")
				
            end

			if flux_scheme == "Roe_FDS"
				roe_fds!(f,rhoL, uL, vL, eL, rhoR, uR, vR, eR,
					jxSA, jySA, SA, view(arr_Fj,:,j,i),
                    QL,QR,FL,FR,LAM,R,Rinv)


			
			else 
				error("ERROR Fluid2d::rhs > Flux scheme cannnot be specified.")				
            end
		end
	end




    @inbounds for i=1:(Ni-2*Nb)#(int i = 0; i < Ni-2*Nb; ++i) {
		for j=1:(Nj-2*Nb) #(int j = 0; j < Nj-2*Nb; ++j) {
			for k=1:Nf #(int k = 0; k < Nf; ++k) {
				arr_Q[k,j,i]= arr_Fi[k,j,i] - arr_Fi[k,j,i+1] + 
							arr_Fj[k,j,i] - arr_Fj[k,j+1,i]; 

            end
		end
	end

end

function reconstruction_basic_muscl3(
    rho0,rho1,
    rho2,rho3,
    u0,u1,
    u2,u3,
    v0,v1,
    v2,v3,
    e0,e1,
    e2,e3
    ) 
    rhoL, rhoR = Fnd.MUSCL3(rho0, rho1, rho2, rho3)
    uL, uR = Fnd.MUSCL3(u0, u1, u2, u3)
	vL, vR = Fnd.MUSCL3(v0, v1, v2, v3)
	eL, eR = Fnd.MUSCL3(e0, e1, e2, e3)
    return rhoL, rhoR,uL, uR,vL, vR ,eL, eR 
end

function reconstruction_basic_mp5(
    rho0,rho1,
    rho2,rho3,
    rho4,rho5,
    u0,u1,
    u2,u3,
    u4,u5,
    v0,v1,
    v2,v3,
    v4,v5,
    e0,e1,
    e2,e3,
    e4,e5
    ) 
    rhoL, rhoR = MP5(rho0, rho1, rho2, rho3, rho4, rho5 )
	uL, uR = MP5(u0, u1, u2, u3, u4, u5)
	vL, vR = MP5(v0, v1, v2, v3, v4, v5)
	eL, eR = MP5(e0, e1, e2, e3, e4, e5)

    return rhoL,uL,vL,eL,rhoR,uR,vR,eR
end