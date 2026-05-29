using KrakenFortran
using TestItems

@testitem "Core kraken calculation" tags = [:unit, :core] begin
    frq = 15.0
    nm = 5

    b1 = [2000 0.0 5000.000]
    ssp1 = [
        0.000 1509.520 0.000 1.03000 0.000 0.000
        10.000 1509.500 0.000 1.03000 0.000 0.000
        20.000 1509.490 0.000 1.03000 0.000 0.000
        5000.000 1536.150 0.000 1.03000 0.000 0.000
    ]

    b2 = [300.0 0.0 5300.0]
    ssp2 = [
        5000.000 1513.108 0.000 1.75100 0.000 0.000
        5300.000 1613.108 0.000 1.75100 0.000 0.000
    ]

    zrc = [0.0, 5300.0]
    zsr = 500.0

    b = [b1; b2]
    ssp = [ssp1; ssp2]

    sspTHS = [0.0 1700.0 0.0 1.7 0.0 0.0]
    sspBHS = [5300.0 1749.0 0.0 1.941 0.1 0.0]
    sspHS = [sspTHS; sspBHS]

    env = UnderwaterEnvironmentFORTRAN(ssp, b, sspHS, zrc, zsr)

    res = kraken(env, frq; n_modes = nm)

    @test haskey(res, "cg")
    @test haskey(res, "cp")
    @test haskey(res, "kr_real")
    @test haskey(res, "kr_imag")
    @test haskey(res, "zm")
    @test haskey(res, "modes")

    @test size(res["modes"], 2) == nm
    @test !all(iszero, res["modes"])
end
