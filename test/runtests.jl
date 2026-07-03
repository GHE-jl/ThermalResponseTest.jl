using ThermalResponseTest
using GroundHeatExchanger: ils, convolution
using Test

@testset "ThermalResponseTest.jl" begin
    # Synthetic infinite-line-source TRT with known parameters (uniform time steps).
    rb, Cs = 0.075, 2.0e6
    k_true, Rb_true, T0 = 2.6, 0.11, 10.0
    t = collect(60.0:60.0:7 * 24 * 3600.0)          # 1 week, 1-min steps
    q = fill(45.0, length(t))                        # constant 45 W/m
    g = ils(t, rb, k_true, Cs)
    T = T0 .+ Rb_true .* q .+ convolution(q, g)

    @testset "FOA — heating temperature (UFOA-T-H)" begin
        k, Rb, reg, ind = fit_ils_foa_T(t, T, q, rb, Cs)
        @test isapprox(k, k_true; rtol = 0.05)
        @test isapprox(Rb, Rb_true; atol = 0.02)
        @test size(reg, 2) == 2
    end

    @testset "FOA — recovery temperature (UFOA-T-R)" begin
        # Build a heating + recovery signal by superposition (heat off after t̄).
        t̄ = 3 * 24 * 3600.0
        qHR = [tᵢ <= t̄ ? 45.0 : 0.0 for tᵢ in t]
        THR = T0 .+ qHR .* Rb_true .+ convolution(qHR, ils(t, rb, k_true, Cs))
        rec = t .> t̄
        k, T0est, reg, ind = fit_ils_foa_T_recovery(t[rec], THR[rec], 45.0, t̄)
        @test isapprox(k, k_true; rtol = 0.1)
    end

    @testset "Model inversion — ILS recovers k and Rb" begin
        res = fit_ils(t, T, q, rb, T0, Cs)
        @test isapprox(res.k, k_true; rtol = 0.02)
        @test isapprox(res.Rb, Rb_true; atol = 0.01)
    end

    @testset "utilities" begin
        @test mean_fluid_temperature([10.0], [12.0], :arithmetic)[1] ≈ 11.0
        @test critical_time(0.075, 2.5, 2.0e6) ≈ 5 * 0.075^2 / (2.5 / 2.0e6)
        @test residence_time(1.0, 100.0) ≈ 200.0
    end
end
