# using Parameters
using UnPack: UnPack, @unpack

# --- Library Handling ---
# During development, we use the local library.
# Once KrakenFortran_jll is registered, this will be replaced by:
# using KrakenFortran_jll
# const libkraken = KrakenFortran_jll.libkraken

function get_libpath()
    lib_ext = Sys.isapple() ? "dylib" : Sys.iswindows() ? "dll" : "so"
    return joinpath(@__DIR__, "kraken.$lib_ext")
end

const libkraken = get_libpath()
# ------------------------

fillnan(x) = isnan(x) ? zero(x) : x

###########################################################################################
# ABSTRACT TYPES
###########################################################################################


###########################################################################################
# STRUCTS
###########################################################################################

"""
`env = Env(ssp, layers, sspHS, z_krak, z_source; n_layers=2, note1="NVW", note2="A", bsig=0, z_step=1.0)`

Creates an underwater environment in the form of a sound speed profile (`ssp`), a
half-space sound speed profile (`sspHS`), and a bottom profile (`bottom`).

The environment is used to compute the propagation of sound waves using Kraken.

Required arguments:
- `ssp`: sound speed profile
- `layers`: bottom profile
- `sspHS`: half-space sound speed profile
- `z_krak`: depth of receivers
- `z_source`: depth of source


Optional keywords:
- `n_layers`: number of layers in the environment (default: 2)
- `note1`: note for the sound speed profile (default: "NVW")
- `note2`: note for the half-space sound speed profile (default: "A")
- `bsig`: bottom interfacial roughness (default: 0)
- `n_sr`: number of sources (default: 1)
- `z_sr`: depth of sources (default: 500)
- `n_bc`: number of boundary conditions (default: size(ssp, 1))
"""
struct UnderwaterEnvironmentFORTRAN
    n_layers::Int
    note1::String
    note2::String
    bsig::Int# bottom interfacial roughness

    ssp::Matrix{Float64}
    sspHS::Matrix{Float64}
    b::Matrix{Float64}

    n_source::Int
    z_source::Float64

    n_krak::Int
    z_krak::Array{Float64,1}
    n_bc::Int
end

function UnderwaterEnvironmentFORTRAN(
    ssp,
    b,
    sspHS,
    z_krak,
    z_source;
    note1 = "NVW",
    note2 = "A",
    bsig = 0,
    z_step = 1.0,
)
    n_layers = size(b, 1)
    n_source = length(z_source)
    n_krak = length(range(z_krak[1], z_krak[2]; step = z_step))
    n_bc = size(ssp, 1)

    if ssp[1, 4] > 500.025  # if ever I use an ssp that conforms to SI units
        ssp[:, 4] /= 1000.0
        sspHS[2, 4] /= 1000.0
    end
    return UnderwaterEnvironmentFORTRAN(
        n_layers,
        note1,
        note2,
        bsig,
        ssp,
        sspHS,
        b,
        n_source,
        z_source,
        n_krak,
        z_krak,
        n_bc,
    )
end

"""
            n_zr = 1
        ρ1 = 1.6, α1 = 0.025, cb = 1900.0, ρb = 2.0, αb = 0.25, type = "1layer_constant")

Helper function to build an underwater environment for use with creating an Env object.

Optional keywords:
- `hw`: water depth (default: 71)
- `cw`: water sound speed (default: 1471)
- `ρw`: water density (default: 1)
- `αw`: water attenuation (default: 0)
- `h1`: sediment depth (default: 10)
- `c1`: sediment sound speed (default: 1500)
- `ρ1`: sediment density (default: 1.6)
- `α1`: sediment attenuation (default: 0.025)
- `cb`: bottom sound speed (default: 1900)
- `ρb`: bottom density (default: 2)
- `αb`: bottom attenuation (default: 0.25)
- `type`: type of environment to build (default: "1layer_constant")

    - "1layer_constant": 1 layer, constant sound speed
"""
function env_builder(;
    hw = 71.0,
    cw = 1471.0,
    ρw = 1.0,
    αw = 0.0,
    h1 = 10.0,
    c1 = 1500.0,
    ρ1 = 1.6,
    α1 = 0.05,
    cb = 1900.0,
    ρb = 2.0,
    αb = 0.25,
    type = "1layer_constant",
)
    if type == "1layer_constant"
        @assert length(cw) == 1
        d0 = hw
        d1 = hw + h1

        b0 = [0.0 0.0 d0]
        ssp0 = [
            0.0 cw 0.0 ρw αw 0.0
            d0 cw 0.0 ρw αw 0.0
        ]

        b1 = [0.0 0.0 d1]
        ssp1 = [
            d0 c1 0.0 ρ1 α1 0.0
            d1 c1 0.0 ρ1 α1 0.0
        ]

        ssp_bhs = [d1 cb 0.0 ρb αb 0.0]
        ssp_ths = [0.0 343.0 0.0 0.00121 0.0 0.0]

        ssp = [ssp0; ssp1]
        sspHS = [ssp_ths; ssp_bhs]
        b = [b0; b1]
    elseif type == "1layer_constant_variable_ssp"
        @assert length(cw) == 4
        d0 = hw
        d1 = hw + h1

        b0 = [0.0 0.0 d0]
        ssp0 = [
            0.0 cw[1] 0.0 ρw αw 0.0
            20.0 cw[2] 0.0 ρw αw 0.0
            40.0 cw[3] 0.0 ρw αw 0.0
            d0 cw[4] 0.0 ρw αw 0.0
        ]

        b1 = [0.0 0.0 d1]
        ssp1 = [
            d0 c1 0.0 ρ1 α1 0.0
            d1 c1 0.0 ρ1 α1 0.0
        ]

        ssp_bhs = [d1 cb 0.0 ρb αb 0.0]
        ssp_ths = [0.0 343.0 0.0 0.00121 0.0 0.0]

        ssp = [ssp0; ssp1]
        sspHS = [ssp_ths; ssp_bhs]
        b = [b0; b1]
    else
        error("type not recognized")
    end

    return ssp, b, sspHS
end

"""
    kraken(env::Env, freq=15.0; n_modes=5, range_max=5e3, c_low=0.0, c_high=nothing)

Call KRAKEN to compute the propagation of sound waves in an underwater environment.

Required arguments:
- `env`: underwater environment

Optional keywords:
- `freq`: frequency (default: 15)
- `n_modes`: number of modes (default: 5)
- `range_max`: maximum range (default: 5e3)
- `c_low`: minimum sound speed (default: 0)
- `c_high`: maximum sound speed (default: maximum of SSP and SSPHS)
"""
function kraken(
    env::UnderwaterEnvironmentFORTRAN,
    freq = 15.0;
    n_modes = 5,
    range_max = 1e4,
    c_low = 0.0,
    c_high = nothing,
)
    if c_high === nothing
        c_high = maximum([maximum(env.ssp[:, 2]), maximum(env.sspHS[:, 2])])
    end

    c_low_high = [c_low c_high]
    @unpack n_layers,
    note1,
    n_bc,
    note2,
    bsig,
    ssp,
    sspHS,
    b,
    n_source,
    z_source,
    n_krak,
    z_krak = env

    cg, cp, kr_real, kr_imag, zm, modes = call_kraken(
        n_modes,
        freq,
        n_layers,
        [UInt8(x) for x in env.note1],
        b,
        n_bc,
        ssp,
        [UInt8(x) for x in env.note2],
        bsig,
        sspHS,
        c_low_high,
        range_max,
        n_source,
        z_source,
        n_krak,
        z_krak,
    )
    return Dict(
        "cg" => cg,
        "cp" => cp,
        "kr_real" => kr_real,
        "kr_imag" => kr_imag,
        "zm" => zm,
        "modes" => modes,
    )
end

"""
    call_kraken(nm, frq, nl, note1, b, nc, ssp, note2, bsig, sspHS, clh, rng, nsr, zsr, nrc, zrc)

Direct call to KRAKEN library.
!!! To be used with `kraken` function only. Not recommened to use as is.

Required arguments:
- `nm`: number of modes
- `frq`: frequency
- `nl`: number of layers
- `note1`: note for SSP
- `b`: bottom profile
- `nc`: number of boundary conditions
- `ssp`: SSP
- `note2`: note for SSPHS
- `bsig`: bottom interfacial roughness
- `sspHS`: SSPHS
- `clh`: minimum and maximum sound speed
- `rng`: maximum range
- `nsr`: number of sources
- `zsr`: depth of sources
- `nrc`: number of receivers
- `zrc`: depth of receivers

"""
function call_kraken(
    nm,
    frq,
    nl,
    note1,
    b,
    nc,
    ssp,
    note2,
    bsig,
    sspHS,
    clh,
    rng,
    nsr,
    zsr,
    nrc,
    zrc,
)
    nz = nsr + nrc

    cg = zeros(1, nm)
    cp = zeros(1, nm)
    kr_real = zeros(1, nm)
    kr_imag = zeros(1, nm)
    zm = zeros(nz, 1)
    modes = zeros(nz, nm)

    ccall(
        (:kraken_, libkraken),
        Nothing,
        (
            Ref{Int},
            Ref{Float64},
            Ref{Int},
            Ref{UInt8},
            Ref{Float64},
            Ref{Int},
            Ref{Float64},
            Ref{UInt8},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Int},
            Ref{Float64},
            Ref{Int},
            Ref{Float64},
            Ref{Int},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
        ),
        nm,
        frq,
        nl,
        note1,
        b,
        nc,
        ssp,
        note2,
        bsig,
        sspHS,
        clh,
        rng,
        nsr,
        zsr,
        nrc,
        zrc,
        nz,
        cg,
        cp,
        kr_real,
        kr_imag,
        zm,
        modes,
    )
    return cg, cp, kr_real, kr_imag, zm, modes
end
