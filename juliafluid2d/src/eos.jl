function calc_eigen(f::Fluid2d{Ni,Nj,Nb,Nf},
    rho,u,v,
    e,ix,iy,
    LAM,R,Rinv
    ) where {Ni,Nj,Nb,Nf}
    error("$(typeof(f)) is not supported")
end

function calc_eigen!(f::IdealGas2d{Ni,Nj,Nb,Nf},
    rho,u,v,
    e,ix,iy,
    LAM,R,Rinv
    ) where {Ni,Nj,Nb,Nf}

    cs = calc_cs(f,rho, u, v, e)
    p = calc_p(f,rho, u, v, e)
    h = (e + p) / rho
    sqr = sqrt(ix * ix + iy * iy)
    ixb = ix / sqr
	iyb = iy / sqr
	bigU = ix * u + iy * v
	bigUb = bigU / sqr
	b1 = 0.5 * (u * u + v * v) * (f.gamma - 1.0) / cs / cs
	b2 = (f.gamma - 1.0) / cs / cs

    LAM[1] = bigU - cs * sqr
	LAM[2] = bigU
	LAM[3] = bigU + cs * sqr
	LAM[4] = bigU

    R[0+1,0+1] = 1.0
	R[1+1,0+1] = 1.0
	R[2+1,0+1] = 1.0
	R[3+1,0+1] = 0.0
	R[0+1,1+1] = u - ixb * cs
	R[1+1,1+1] = u
	R[2+1,1+1] = u + ixb * cs
	R[3+1,1+1] = -iyb
	R[0+1,2+1] = v - iyb * cs
	R[1+1,2+1] = v
	R[2+1,2+1] = v + iyb * cs
	R[3+1,2+1] = ixb
	R[0+1,3+1] = h - cs * bigUb
	R[1+1,3+1] = 0.5 * (u * u + v * v)
	R[2+1,3+1] = h + cs * bigUb
	R[3+1,3+1] = -(iyb * u - ixb * v)

    Rinv[0+1,0+1] = 0.5 * (b1 + bigUb / cs)
	Rinv[1+1,0+1] = -0.5 * (ixb / cs + b2 * u)
	Rinv[2+1,0+1] = -0.5 * (iyb / cs + b2 * v)
	Rinv[3+1,0+1] = 0.5 * b2
	Rinv[0+1,1+1] = 1.0 - b1
	Rinv[1+1,1+1] = b2 * u
	Rinv[2+1,1+1] = b2 * v
	Rinv[3+1,1+1] = -b2
	Rinv[0+1,2+1] = 0.5 * (b1 - bigUb / cs)
	Rinv[1+1,2+1] = 0.5 * (ixb / cs - b2 * u)
	Rinv[2+1,2+1] = 0.5 * (iyb / cs - b2 * v)
	Rinv[3+1,2+1] = 0.5 * b2
	Rinv[0+1,3+1] = iyb * u - ixb * v
	Rinv[1+1,3+1] = -iyb
	Rinv[2+1,3+1] = ixb
	Rinv[3+1,3+1] = 0.0




end