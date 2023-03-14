
include("fluid2d.jl")

using .Fluid2d_module

function main()
    Ni = 108
    Nj = 108
    Nb = 4
    Nf = 4
    Tmax = 3
    Nout = 100

    dir = "../output/"
    dir_o = dir*"data_grid108x108/"
    f_coordinate = dir*"coordinate_grid108x108.dat"
    f_settings = dir*"settings.dat"

    fluid = IdealGas2d(Ni,Nj,Nb,Nf,dir_o,f_coordinate,f_settings)
    iter = 0;
	t = 0.0;
	Tout = Tmax / Nout;

    @time for tstep=1:Nout
        iter = 0
       while t < tstep *Tout
            dt = calc_cfl(fluid,0.7)
            march_ssprk3!(fluid,dt, "periodical_in_i", "MP5_basic", "Roe_FDS")
            t += dt
            iter += 1
            #println("$t $iter $dt")
        end
        output_basic(fluid,tstep, 
            iter
            ) 
    end
end

main()