function calc_conserved!(rho,u,v,e,S,Q)
    Q[1] = rho * S
	Q[2] = rho * u * S
	Q[3] = rho * v * S
	Q[4] = e * S
end

function roe_fds!(f::Fluid2d{Ni,Nj,Nb,Nf},
    rhoL,uL,vL,
    eL,rhoR,uR,
    vR,eR,ixS,
    iyS,S,Fc,
    QL,QR,FL,FR,LAM,R,Rinv
    ) where {Ni,Nj,Nb,Nf}
    ix = ixS / S
    iy = iyS / S
    QL .= 0
    QR .= 0
    FL .= 0
    FR.=0
    LAM .= 0
    R .= 0
    Rinv .=0
    rhoA, uA, vA, eA = roe_average(f,rhoL, uL, vL, eL, rhoR, uR, vR, eR);
    
    #QL = zeros(Nf)
    #QR = zeros(Nf)
    calc_conserved!(rhoL, uL, vL, eL, S, QL)
    calc_conserved!(rhoR, uR, vR, eR, S, QR)
    #FL = zeros(Nf)
    #FR = zeros(Nf)
    calc_flux_conv!(f,rhoL, uL, vL, eL, ixS, iyS, FL);
	calc_flux_conv!(f,rhoR, uR, vR, eR, ixS, iyS, FR);

    #LAM = zeros(Nf)
    #R = zeros(Nf,Nf)
    #Rinv = zeros(Nf,Nf)
	calc_eigen!(f,rhoA, uA, vA, eA, ix, iy, LAM, R, Rinv)
    #error()

    @inbounds for i=1:Nf#(int i = 0; i < Nf; ++i) {
		Fc[i] = FR[i] + FL[i];
		for j=1:Nf#(int j = 0; j < Nf; ++j) {
			for k=1:Nf #(int k = 0; k < Nf; ++k) {
				Fc[i] -= R[j,i] * abs(LAM[j]) * Rinv[k,j] * (QR[k] - QL[k])
            end
		end
		Fc[i] *= 0.5;
	end

end

function roe_average(f::Fluid2d{Ni,Nj,Nb,Nf},
    rhoL,uL,vL,
    eL,rhoR,uR,
    vR,eR
    ) where {Ni,Nj,Nb,Nf}

    pL = calc_p(f,rhoL, uL, vL, eL)
    pR = calc_p(f,rhoR, uR, vR, eR)
    hL = (eL + pL) / rhoL
    hR = (eR + pR) / rhoR
    sqrL = sqrt(rhoL)
	sqrR = sqrt(rhoR)

    rhoA = sqrL * sqrR
	uA = (sqrL * uL + sqrR * uR) / (sqrL + sqrR)
	vA = (sqrL * vL + sqrR * vR) / (sqrL + sqrR)

    hA = (sqrL * hL + sqrR * hR) / (sqrL + sqrR)
    eA = calc_e(f,rhoA, uA, vA, hA)

    return rhoA,uA,vA,eA
end

function calc_flux_conv!(f::Fluid2d{Ni,Nj,Nb,Nf},
    rho,u,
    v,e,ixS,iyS,Fc
    ) where {Ni,Nj,Nb,Nf}
    bigU = ixS * u + iyS * v
    p = calc_p(f,rho, u, v, e)

    Fc[1] = rho * bigU
	Fc[2] = rho * u * bigU + ixS * p
	Fc[3] = rho * v * bigU + iyS * p
	Fc[4] = (e + p) * bigU
end

function calc_basic(f::Fluid2d{Ni,Nj,Nb,Nf},Q,S) where {Ni,Nj,Nb,Nf}
        rho = Q[1] / S;
        u = Q[2] / rho / S;
        v = Q[3] / rho / S;
        e = Q[4] / S;
        return rho,u,v,e
end