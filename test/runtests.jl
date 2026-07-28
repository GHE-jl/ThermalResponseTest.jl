using ThermalResponseTest
using GroundHeatExchanger: ils, convolution, ILSModel
using DataFrames
using Test

@testset "ThermalResponseTest.jl" begin
    # Synthetic infinite-line-source TRT with known parameters (uniform time steps).
    rb, Cs = 0.075, 2.0e6
    k_true, Rbₑ_true, T0 = 2.6, 0.11, 10.0
    t = collect(60.0:60.0:7 * 24 * 3600.0)          # 1 week, 1-min steps
    q = fill(45.0, length(t))                        # constant 45 W/m
    g = ils(t, rb, k_true, Cs)
    T = T0 .+ Rbₑ_true .* q .+ convolution(q, g)

    @testset "FOA — heating temperature (UFOA-T-H)" begin
        k, Rbₑ, reg, ind = fit_ils_foa_T(t, T, q, rb, T0, Cs)
        @test isapprox(k, k_true; rtol = 0.05)
        @test isapprox(Rbₑ, Rbₑ_true; atol = 0.02)
        @test size(reg, 2) == 2
    end

    @testset "FOA — recovery temperature (UFOA-T-R)" begin
        # Build a heating + recovery signal by superposition (heat off after t̄).
        t̄ = 3 * 24 * 3600.0
        qHR = [tᵢ <= t̄ ? 45.0 : 0.0 for tᵢ in t]
        THR = T0 .+ qHR .* Rbₑ_true .+ convolution(qHR, ils(t, rb, k_true, Cs))
        rec = t .> t̄
        k, reg, ind = fit_ils_foa_T_recovery(t[rec], THR[rec], 45.0, rb, t̄, Cs)
        @test isapprox(k, k_true; rtol = 0.1)
    end

    @testset "Model inversion — ILS recovers k and Rbₑ" begin
        res = fit_ils(t, T, q, rb, T0, Cs)
        @test isapprox(res.k, k_true; rtol = 0.02)
        @test isapprox(res.Rbₑ, Rbₑ_true; atol = 0.01)
    end

    @testset "fit_ground_response — equivalence with the fit_ils wrapper" begin
        # fit_ils is a thin wrapper around fit_ground_response with an ILSModel; both must agree.
        res_direct = fit_ils(t, T, q, rb, T0, Cs)
        res_generic = fit_ground_response(t, T, q, rb, T0, p -> ILSModel(p[1], Cs),
            [2.5], [0.2], [7.0])
        @test isapprox(res_direct.k, res_generic.params[1]; rtol = 1e-6)
        @test isapprox(res_direct.Rbₑ, res_generic.Rbₑ; rtol = 1e-6)
    end

    @testset "decompose_trt / TRTDataset — recirculation phase before heating" begin
        # A recirculation phase (idle pump, power = 0, T ≈ T0) logged before heating starts must not
        # bias any TRTDataset-based fit: :t_rel rebases time to t = 0 at the start of heating.
        Δt = 600.0                                        # 10-minute logging step
        t_local = collect(Δt:Δt:3 * 24 * 3600.0)          # 3 days of heating
        t̄ = 1 * 24 * 3600.0                               # heating stops after 1 day
        q_local = [tᵢ <= t̄ ? 1000.0 : 0.0 for tᵢ in t_local]  # W (H = 1 m, so also W/m)
        g_local = ils(t_local, rb, k_true, Cs)
        T_local = T0 .+ q_local .* Rbₑ_true .+ convolution(q_local, g_local)

        n_recirc = 3                                      # 30 minutes of recirculation
        recirc_time = collect(0.0:Δt:(n_recirc - 1) * Δt)
        trt = DataFrame(
            elapsed_time = [recirc_time; t_local .+ n_recirc * Δt],
            power        = [zeros(n_recirc); q_local],
            T_mean       = [fill(T0, n_recirc); T_local],
        )

        dataset = decompose_trt(trt)
        @test nrow(dataset.heating) == count(>(0), q_local)
        @test nrow(dataset.recovery) == count(==(0), q_local)
        @test dataset.heating.t_rel ≈ t_local[q_local .> 0]        # rebased clock == the true one
        @test dataset.recovery.t_rel ≈ t_local[q_local .== 0]

        # Inversion and FOA must recover the true parameters despite the recirculation phase.
        res = fit_ils(dataset, 1.0, rb, T0, Cs)
        @test isapprox(res.k, k_true; rtol = 0.05)
        @test isapprox(res.Rbₑ, Rbₑ_true; atol = 0.02)

        k_H, Rbₑ_H, _, _ = fit_ils_foa_T(dataset, 1.0, rb, Cs)
        @test isapprox(k_H, k_true; rtol = 0.1)
        @test isapprox(Rbₑ_H, Rbₑ_true; atol = 0.02)

        k_R, _, _ = fit_ils_foa_T_recovery(dataset, 1.0, rb, Cs)
        @test isapprox(k_R, k_true; rtol = 0.15)
    end

    @testset "utilities" begin
        @test mean_fluid_temperature([10.0], [12.0], :arithmetic)[1] ≈ 11.0
        @test critical_time(0.075, 2.5, 2.0e6) ≈ 5 * 0.075^2 / (2.5 / 2.0e6)
        @test residence_time(1.0, 100.0) ≈ 200.0
    end
end
