using KrakenFortran
using TestItems

@testitem "Environment builder" tags = [:unit, :env] begin
    ssp, b, sspHS = env_builder()

    @test size(ssp, 1) == 4
    @test size(b, 1) == 2
    @test size(sspHS, 1) == 2

    # Test with variable SSP
    cw_var = [1500.0, 1490.0, 1480.0, 1470.0]
    ssp_v, b_v, sspHS_v = env_builder(cw = cw_var, type = "1layer_constant_variable_ssp")
    @test size(ssp_v, 1) == 6
end

@testitem "UnderwaterEnvironmentFORTRAN constructor" tags = [:unit, :env] begin
    ssp, b, sspHS = env_builder()
    zrc = [0.0, 100.0]
    zsr = 50.0

    env = UnderwaterEnvironmentFORTRAN(ssp, b, sspHS, zrc, zsr)

    @test env.n_layers == 2
    @test env.n_source == 1
    @test env.z_source == 50.0
end
