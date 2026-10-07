using Test
using DataFrames: nrow
using Statistics: mean
using SyntheticReefDomain
using SyntheticReefDomain: attr_spec, resolve_attr, climate_params, TOPOLOGIES, LAYOUTS
using ADRIA

# Set RUN_ADRIA_SCENARIOS=false to skip the (slow) ADRIA scenario runs.
const RUN_SCENARIOS = get(ENV, "RUN_ADRIA_SCENARIOS", "true") == "true"

const TMP = mktempdir()
ADRIA.setup()  # sets ADRIA_THRESHOLD etc., needed by run_scenarios
ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

small_cfg(; n=10, name="Test$(n)", kwargs...) = DomainConfig(;
    name=name, n_dhw_scenarios=6, n_cyclone_scenarios=5,
    reefs=ReefConfig(; n=n), kwargs...
)

row_sums(C) = vec(sum(C; dims=2))

@testset "SyntheticReefDomain" begin
    @testset "Configuration" begin
        cfg = DomainConfig()
        @test validate(cfg) === cfg

        # TOML round trip
        path = joinpath(TMP, "roundtrip.toml")
        save_config(path, cfg)
        @test config_to_dict(load_config(path)) == config_to_dict(cfg)

        # The shipped configs load and validate
        cfg_dir = joinpath(@__DIR__, "..", "configs")
        for fn in [joinpath(cfg_dir, "default.toml");
                   readdir(joinpath(cfg_dir, "examples"); join=true)]
            @test validate(load_config(fn)) isa DomainConfig
        end

        # Attribute specifications
        @test attr_spec(3) === 3.0
        @test attr_spec([1, 2]) == [1.0, 2.0]
        @test attr_spec(Dict("min" => 1, "max" => 2)) == (min=1.0, max=2.0)
        @test attr_spec((from=0.5, to=1.5)) == (from=0.5, to=1.5)
        @test_throws ArgumentError attr_spec("high")
        @test resolve_attr((from=0.0, to=1.0), 3, nothing) == [0.0, 0.5, 1.0]
        @test resolve_attr((from=0.0, to=1.0), 1, nothing) == [0.5]

        # Invalid configurations
        @test_throws ArgumentError validate(DomainConfig(; reefs=ReefConfig(; n=0)))
        @test_throws ArgumentError validate(DomainConfig(; reefs=ReefConfig(; n=21)))
        @test_throws ArgumentError validate(DomainConfig(; timeframe=[2025, 2040]))
        @test_throws ArgumentError validate(DomainConfig(; name="bad name"))
        @test_throws ArgumentError validate(
            DomainConfig(; reefs=ReefConfig(; n=3, area_m2=[1e4, 2e4]))
        )
        @test_throws ArgumentError validate(
            DomainConfig(; connectivity=ConnectivityConfig(; topology="star"))
        )
        @test_throws ArgumentError load_config(Dict("reefs" => Dict("nn" => 3)))
        @test_throws ArgumentError load_config(
            Dict("dhw" => Dict("climate" => Dict("045" => Dict())))
        )
        @test_throws ArgumentError load_config(
            Dict("dhw" => Dict("climate" => Dict("low" => Dict())))
        )

        # Climate overrides
        dc = DHWConfig(; climate=Dict(1 => Dict("trend_per_decade" => 0.5), 2 => Dict()))
        @test climate_params(dc, 1).trend_per_decade == 0.5
        @test climate_params(dc, 2).trend_per_decade == dc.trend_per_decade
    end

    @testset "Layout and reef attributes" begin
        for layout in LAYOUTS, n in (1, 2, 7, 20)
            cfg = small_cfg(; n=n, reefs=ReefConfig(; n=n, layout=layout, n_clusters=min(n, 3)))
            reefs = build_reefs(cfg)
            @test nrow(reefs) == n
            @test reefs.reef_siteid == [SyntheticReefDomain.reef_id(i) for i in 1:n]
            @test issorted(reefs.reef_siteid)
            @test length(unique(reefs.cluster_id)) == min(n, 3)
            @test all(45 .<= reefs.k .<= 63)
        end

        # Identical (ub_med, depth_med) pairs are nudged apart for C~scape
        cfg = small_cfg(; reefs=ReefConfig(; n=5, ub_med=0.3, depth_med=5.0))
        reefs = build_reefs(cfg)
        @test length(unique(zip(reefs.ub_med, reefs.depth_med))) == 5
        @test all(abs.(reefs.ub_med .- 0.3) .< 1e-3)

        # Adding reefs does not change existing reefs' attributes
        a = build_reefs(small_cfg(; reefs=ReefConfig(; n=5)))
        b = build_reefs(small_cfg(; reefs=ReefConfig(; n=8)))
        @test a.area == b.area[1:5]
        @test a.dhw_exposure == b.dhw_exposure[1:5]
    end

    @testset "Connectivity" begin
        n = 6
        conn(; kw...) = begin
            cfg = small_cfg(; n=n, connectivity=ConnectivityConfig(; kw...))
            build_connectivity(cfg, build_reefs(cfg))
        end
        offdiag(C) = [C[i, j] for i in 1:n, j in 1:n if i != j]

        for topology in TOPOLOGIES
            C = conn(; topology=topology, strength=0.3, self_recruitment=0.25)
            @test size(C) == (n, n)
            @test all(0 .<= C .<= 1)
            @test all(row_sums(C) .<= 1 + 1e-12)
            @test all(C[i, i] == 0.25 for i in 1:n)
        end

        @test all(offdiag(conn(; topology="isolated")) .== 0)
        @test all(offdiag(conn(; topology="all_to_all")) .> 0)

        fwd = conn(; topology="chain", asymmetry=1.0)
        @test all(fwd[i, i + 1] > 0 for i in 1:(n - 1))
        @test all(fwd[i + 1, i] == 0 for i in 1:(n - 1))
        @test count(>(0), fwd) == n + (n - 1)

        sym = conn(; topology="chain", asymmetry=0.0)
        @test sym == sym'

        ring = conn(; topology="ring", asymmetry=1.0)
        @test ring[n, 1] > 0

        reach2 = conn(; topology="chain", asymmetry=1.0, reach=2)
        @test reach2[1, 3] ≈ reach2[1, 2] / 2

        @test all(offdiag(conn(; topology="random", density=0.0)) .== 0)

        # Export budget: self + exports never exceeds 1
        full = conn(; topology="all_to_all", strength=1.0, self_recruitment=0.4)
        @test all(row_sums(full) .≈ 1.0)

        # Distance decay weakens with distance
        dd = conn(; topology="distance_decay", decay_km=5.0)
        @test dd[1, 2] > dd[1, 3] > dd[1, 4]

        # One reef: self recruitment only
        cfg1 = small_cfg(; n=1)
        @test build_connectivity(cfg1, build_reefs(cfg1)) == fill(0.2, 1, 1)
    end

    @testset "Initial cover" begin
        cfg = small_cfg(; cover=CoverConfig(; total_percent=30.0))
        cover = build_cover(cfg, build_reefs(cfg))
        @test size(cover) == (5, 10)
        @test all(vec(sum(cover; dims=1)) .≈ 0.30)

        # Percent of reef area -> fraction of habitable area
        cfg = small_cfg(;
            reefs=ReefConfig(; n=4, k_percent=50.0),
            cover=CoverConfig(; total_percent=20.0, basis="reef")
        )
        @test all(vec(sum(build_cover(cfg, build_reefs(cfg)); dims=1)) .≈ 0.40)

        cfg = small_cfg(;
            reefs=ReefConfig(; n=4, k_percent=50.0),
            cover=CoverConfig(; total_percent=60.0, basis="reef")
        )
        @test_throws ArgumentError build_cover(cfg, build_reefs(cfg))

        # Composition jitter keeps totals and zero-weight groups
        cfg = small_cfg(;
            cover=CoverConfig(;
                total_percent=25.0, composition=[1, 1, 0, 1, 1], jitter=0.3
            )
        )
        cover = build_cover(cfg, build_reefs(cfg))
        @test all(vec(sum(cover; dims=1)) .≈ 0.25)
        @test all(cover[3, :] .== 0)
        @test length(unique(round.(cover[1, :]; digits=8))) > 1
    end

    @testset "Disturbances" begin
        cfg = small_cfg()
        reefs = build_reefs(cfg)
        dhw, regional = build_dhw(cfg, reefs)
        D = dhw[45]
        @test size(D) == (75, 10, 6)
        @test all(D .>= 0)
        @test size(regional[45]) == (75, 6)
        @test mean(D[66:75, :, :]) > mean(D[1:10, :, :])  # warming trend

        # Identical DHW on every reef
        same = small_cfg(; dhw=DHWConfig(; local_sd=0.0, reef_exposure=1.0, reef_offset=0.0))
        Ds = build_dhw(same, build_reefs(same))[1][45]
        @test all(Ds[:, r, :] == Ds[:, 1, :] for r in 2:10)

        # Climate levels share random draws; only their parameters differ
        two = small_cfg(; dhw=DHWConfig(; climate=Dict(1 => Dict(), 2 => Dict())))
        d2, _ = build_dhw(two, build_reefs(two))
        @test d2[1] == d2[2]
        hot = small_cfg(;
            dhw=DHWConfig(; climate=Dict(1 => Dict(), 2 => Dict("trend_per_decade" => 3.0)))
        )
        dh, _ = build_dhw(hot, build_reefs(hot))
        @test mean(dh[2]) > mean(dh[1])

        # Prescribed heatwave years appear in every scenario
        quiet = small_cfg(; dhw=DHWConfig(; event_prob=0.0))
        forced = small_cfg(;
            dhw=DHWConfig(; event_prob=0.0, event_years=[2050], event_magnitude=[10.0, 10.0])
        )
        _, reg0 = build_dhw(quiet, build_reefs(quiet))
        _, reg = build_dhw(forced, build_reefs(forced))
        Δ = reg[45] .- reg0[45]
        @test all(Δ[26, :] .≈ 10.0)
        @test all(Δ[setdiff(1:75, 26), :] .== 0.0)

        # Reproducible
        @test build_dhw(cfg, reefs)[1][45] == D
        other = small_cfg(; seed=99)
        @test build_dhw(other, build_reefs(other))[1][45] != D

        # Cyclones: system-wide, group specific, within [0, 1]
        M, intensity = build_cyclones(cfg, 10)
        @test size(M) == (75, 10, 5, 5)
        @test all(0 .<= M .<= 1)
        @test all(M[:, r, :, :] == M[:, 1, :, :] for r in 2:10)
        @test any(intensity .> 0)
        hit = findfirst(>(0), intensity)
        @test issorted(M[hit[1], 1, :, hit[2]]; rev=true)  # vulnerability order

        off = small_cfg(; cyclones=CycloneConfig(; enabled=false))
        @test all(build_cyclones(off, 10)[1] .== 0)

        fixed = small_cfg(; cyclones=CycloneConfig(; annual_prob=0.0, event_years=[2030]))
        _, fi = build_cyclones(fixed, 10)
        @test all(fi[6, :] .> 0)
        @test count(>(0), fi) == 5

        # Waves: ub_med repeated, with as many scenarios as DHW
        W = build_waves(cfg, reefs)
        @test size(W) == size(D)
        @test W[17, :, 3] == reefs.ub_med
    end

    @testset "ADRIA round trip" begin
        cases = [(n, "distance_decay") for n in (1, 2, 20)]
        append!(cases, [(10, t) for t in TOPOLOGIES])
        for (n, topology) in cases
            name = "RT_$(n)_$(topology)"
            cfg = small_cfg(; n=n, name=name, connectivity=ConnectivityConfig(; topology=topology))
            synth = build_domain(cfg)
            path = write_domain(synth, joinpath(TMP, name))
            @test isfile(joinpath(path, "spatial", "RT_$(n)_$(topology).gpkg"))

            dom = ADRIA.load_domain(path, "45")
            ids = synth.reefs.reef_siteid
            @test dom.loc_ids == ids
            @test dom.loc_data.reef_siteid == ids
            @test Array(dom.conn) == synth.connectivity
            @test Array(dom.dhw_scens) == synth.dhw[45]
            @test Array(dom.wave_scens) == synth.waves
            @test size(dom.cyclone_mortality_scens) == (75, n, 5, 5)
            @test Array(dom.cyclone_mortality_scens) == synth.cyclone_mortality
            @test dom.loc_data.k ≈ synth.reefs.k ./ 100
            @test dom.loc_data.area ≈ synth.reefs.area
            @test dom.loc_data.depth_med ≈ synth.reefs.depth_med
            @test dom.loc_data.ub_med ≈ synth.reefs.ub_med
            @test dom.loc_data.CB_CALIB_GROUPS == synth.reefs.CB_CALIB_GROUPS

            # ADRIA splits each group into size classes; the sums must match
            init = Array(dom.init_coral_cover)  # (5 groups * 7 sizes) x n
            per_group = dropdims(sum(reshape(init, 7, 5, n); dims=1); dims=1)
            @test per_group ≈ synth.cover

            # Writing again needs force=true; force only replaces data packages
            @test_throws ArgumentError write_domain(synth, path)
            @test write_domain(synth, path; force=true) == path
        end

        # Several climate levels -> one DHW file each
        cfg = small_cfg(; name="TwoClimates", dhw=DHWConfig(; climate=Dict(1 => Dict(), 2 => Dict("trend_per_decade" => 3.0))))
        synth = build_domain(cfg)
        path = write_domain(synth, joinpath(TMP, "TwoClimates"))
        @test isfile(joinpath(path, "DHWs", "dhwRCP1.nc"))
        @test isfile(joinpath(path, "DHWs", "dhwRCP2.nc"))
        @test climate_labels(path) == [1, 2]
        @test Array(ADRIA.load_domain(path, "2").dhw_scens) == synth.dhw[2]

        not_a_domain = mkpath(joinpath(TMP, "NotADomain"))
        touch(joinpath(not_a_domain, "keep.txt"))
        @test_throws ArgumentError write_domain(synth, not_a_domain; force=true)
        @test isfile(joinpath(not_a_domain, "keep.txt"))
    end

    if RUN_SCENARIOS
        @testset "ADRIA scenario runs" begin
            for n in (2, 10, 20)
                cfg = small_cfg(; n=n, name="Run_$n")
                path = write_domain(build_domain(cfg), joinpath(TMP, "Run_$n"))
                dom = ADRIA.load_domain(path, "45")
                scens = ADRIA.sample_cf(dom, 4)
                rs = ADRIA.run_scenarios(dom, scens, "45"; show_progress=false)
                rc = ADRIA.metrics.relative_cover(rs)
                @test size(rc, 2) == n
                @test all(isfinite, rc)
                @test rc[end, :, :] != rc[1, :, :]
            end

            # Two climate levels in one run
            cfg = small_cfg(; name="RunTwoClimates", dhw=DHWConfig(; climate=Dict(1 => Dict(), 2 => Dict("trend_per_decade" => 3.0))))
            path = write_domain(build_domain(cfg), joinpath(TMP, "RunTwoClimates"))
            dom = ADRIA.load_domain(path, "1")
            scens = ADRIA.sample_cf(dom, 4)
            rs = ADRIA.run_scenarios(dom, scens, ["1", "2"]; show_progress=false)
            @test size(ADRIA.metrics.relative_cover(rs), 3) == 8
        end
    end
end
