module Fluid2d_module
    export IdealGas2d,calc_cfl,march_ssprk3!,output_basic

    include("Fnd.jl")

abstract type Fluid2d{Ni,Nj,Nb,Nf} end

mutable struct Fluid2d_core{Ni,Nj,Nb,Nf} <: Fluid2d{Ni,Nj,Nb,Nf}
    t::Float64 
    tstep::Int64 
    x::Matrix{Float64}
    y::Matrix{Float64} 
    rho::Matrix{Float64} 
    u::Matrix{Float64} 
    v::Matrix{Float64}
    e::Matrix{Float64}
    arr_Q0::Array{Float64,3} 
    arr_Q1::Array{Float64,3} 
    arr_Q2::Array{Float64,3} 
    arr_Fi::Array{Float64,3} 
    arr_Fj::Array{Float64,3} 
    dir_o::String 
    f_coordinate::String 
    f_settings::String 
    S::Matrix{Float64} 
    ixS::Matrix{Float64} 
    iyS::Matrix{Float64} 
    jxS::Matrix{Float64} 
    jyS::Matrix{Float64} 
    dx::Matrix{Float64}

    function Fluid2d_core(Ni,Nj,Nb,Nf,
            dir_o,f_coordinate,f_settings)
        f = new{Ni,Nj,Nb,Nf}()
        f.t = 0.0
        f.tstep= 0
        f.x= zeros(Nj,Ni)
        f.y = zeros(Nj,Ni)
        f.rho= zeros(Nj,Ni)
        f.u= zeros(Nj,Ni)
        f.v= zeros(Nj,Ni)
        f.e= zeros(Nj,Ni)
        f.arr_Q0 = zeros(Nf,Nj-2*Nb,Ni-2*Nb)
        f.arr_Q1 = zeros(Nf,Nj-2*Nb,Ni-2*Nb)
        f.arr_Q2= zeros(Nf,Nj-2*Nb,Ni-2*Nb)
        f.arr_Fi= zeros(Nf,Nj-2*Nb,Ni-2*Nb+1)
        f.arr_Fj= zeros(Nf,Nj-2*Nb+1,Ni-2*Nb)
        f.dir_o =  dir_o
        f.f_coordinate = f_coordinate
        f.f_settings = f_settings
        f.S= zeros(Nj-1,Ni-1)
        f.ixS= zeros(Nj-2*Nb+2,Ni-2*Nb+2)
        f.iyS= zeros(Nj-2*Nb+2,Ni-2*Nb+2)
        f.jxS = zeros(Nj-2*Nb+2,Ni-2*Nb+2)
        f.jyS= zeros(Nj-2*Nb+2,Ni-2*Nb+2)
        f.dx= zeros(Nj-2*Nb,Ni-2*Nb)
        initialize!(f,dir_o,f_coordinate,f_settings)
        return f
    end

end

struct IdealGas2d{Ni,Nj,Nb,Nf} <: Fluid2d{Ni,Nj,Nb,Nf}
    fluid::Fluid2d_core{Ni,Nj,Nb,Nf}
    gamma::Float64
    function IdealGas2d(Ni,Nj,Nb,Nf,
        dir_o,f_coordinate,f_settings;gamma = 1.4)
        fluid = Fluid2d_core(Ni,Nj,Nb,Nf,dir_o,f_coordinate,f_settings)
        return new{Ni,Nj,Nb,Nf}(fluid,gamma)
    end
end

function initialize!(f,dir_o,f_coordinate,f_settings)
    f.t = 0
    f.tstep =0
    input_coordinate!(f,f_coordinate)
    input_basic!(f,dir_o*"b0000000.dat")
    calc_metrices_dx!(f)
end

calc_cs(f::Fluid2d,rho,u,v, e) = error("type $(typeof(f)) is not supported!")
calc_cs(f::IdealGas2d,rho,u,v, e) = sqrt(f.gamma * (f.gamma - 1.0) * (e / rho - 0.5 * (u * u + v * v)));

calc_e(f::Fluid2d,rho, u, v, h) = error("type $(typeof(f)) is not supported!")
calc_e(f::IdealGas2d,rho, u, v, h) =rho * (h + (f.gamma - 1.0) * 0.5 * (u * u + v * v)) / f.gamma

calc_p(f::Fluid2d,rho,u,v,e) = error("type $(typeof(f)) is not supported!")
calc_p(f::IdealGas2d,rho,u,v,e) = (f.gamma - 1.0) * (e - 0.5 * rho * (u * u + v * v))


get_rho(f::IdealGas2d) = f.fluid.rho
get_u(f::IdealGas2d)= f.fluid.u
get_v(f::IdealGas2d)= f.fluid.v
get_e(f::IdealGas2d)= f.fluid.e
get_dx(f::IdealGas2d)= f.fluid.dx
get_S(f::IdealGas2d)  = f.fluid.S
get_arr_Q0(f::IdealGas2d)= f.fluid.arr_Q0
get_arr_Q1(f::IdealGas2d)= f.fluid.arr_Q1
get_arr_Q2(f::IdealGas2d)= f.fluid.arr_Q2
get_ixS(f::IdealGas2d)= f.fluid.ixS
get_iyS(f::IdealGas2d)= f.fluid.iyS
get_jxS(f::IdealGas2d)= f.fluid.jxS
get_jyS(f::IdealGas2d)= f.fluid.jyS
get_arr_Fi(f::IdealGas2d)= f.fluid.arr_Fi
get_arr_Fj(f::IdealGas2d)= f.fluid.arr_Fj
get_dir_o(f::IdealGas2d)= f.fluid.dir_o

function add_t!(f::IdealGas2d,dt)
    f.fluid.t += dt
end

function add_tstep!(f::IdealGas2d,dt)
    f.fluid.tstep += dt
end

include("eos.jl")
include("io_settings.jl")
include("marching.jl")
include("flux_scheme.jl")
include("bc.jl")

end



