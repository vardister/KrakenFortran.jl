using KrakenFortran
using TestItems
using NamedArrays

@testitem "Pressure field calculation" tags = [:unit, :field] begin
    ssp, b, sspHS = env_builder()
    zrc = [0.0, 80.0] # water depth is 71 by default + 10 sediment
    zsr = 50.0

    env = UnderwaterEnvironmentFORTRAN(ssp, b, sspHS, zrc, zsr)

    freq = 100.0
    r_vec = [1000.0, 2000.0]
    zr_vec = [10.0, 30.0, 50.0]

    # Complex pressure field
    pf = pressure_field_fortran(env, freq, r_vec, zsr, zr_vec; n_modes = 5)

    @test size(pf) == (3, 2)
    @test eltype(pf) == ComplexF64
    @test !all(iszero, pf)

    # Transmission Loss (TL)
    tl = pressure_field_fortran(env, freq, r_vec, zsr, zr_vec; n_modes = 5, TL = true)

    @test size(tl) == (3, 2)
    @test eltype(tl) == Float64
    @test !all(iszero, tl)
end
