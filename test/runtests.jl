using Test
using DataFrames: nrow
using Statistics: mean
using SyntheticReefDomain
using SyntheticReefDomain: attr_spec, resolve_attr, climate_params, TOPOLOGIES, LAYOUTS
using ADRIA
using CairoMakie  # loads the plotting extension (domain README figures)

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
        @test config_to_dict(load_config(IOBuffer(read(path)))) == config_to_dict(cfg)
        cfg_td = DomainConfig(; connectivity=ConnectivityConfig(;
            topology="truncated_dilution", threshold_m=100.0, dispersal_km=12.5, matrix_total=1.0))
        path_td = joinpath(TMP, "roundtrip_td.toml")
        save_config(path_td, cfg_td)
        @test config_to_dict(load_config(path_td)) == config_to_dict(cfg_td)
        @test load_config(path_td).connectivity.dispersal_km == 12.5

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
        for kw in ((dispersal_km=0.0,), (threshold_m=-1.0,), (matrix_total=0.0,), (matrix_total=1.5,))
            @test_throws ArgumentError validate(DomainConfig(;
                connectivity=ConnectivityConfig(; topology="truncated_dilution", kw...)))
        end
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
            # truncated_dilution does not use self_recruitment (tested below)
            topology == "truncated_dilution" || @test all(C[i, i] == 0.25 for i in 1:n)
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

    @testset "Truncated dilution connectivity" begin
        td(n, reefs=ReefConfig(; n=n); kw...) = begin
            cfg = small_cfg(; n=n, reefs=reefs,
                connectivity=ConnectivityConfig(; topology="truncated_dilution", kw...))
            build_connectivity(cfg, build_reefs(cfg))
        end

        # Default reefs (radius 93-137 m) retain nothing; every pair is linked; scaled total
        n = 6
        C = td(n)
        @test all(C .>= 0)
        @test all(C[i, i] == 0 for i in 1:n)
        @test SyntheticReefDomain.offdiag_links(C) == n * (n - 1)
        @test sum(C) ≈ 0.00262
        C1 = td(n; matrix_total=1.0)
        @test sum(C1) ≈ 1.0
        @test C ≈ 0.00262 .* C1
        @test td(n; dispersal_km=5.0)[1, 3] < C[1, 3]  # shorter dispersal: weaker far links

        # Matches a direct numerical integral of the truncated kernel. Reefs 300 m apart
        # with radii of 169-219 m: local retention on every reef, and the 157.57 m
        # truncation cuts into the links between neighbours.
        reefs3 = ReefConfig(; n=3, spacing_km=0.3, area_m2=[150_000.0, 200_000.0, 250_000.0],
            k_percent=60.0)
        cfg3 = small_cfg(; n=3, reefs=reefs3,
            connectivity=ConnectivityConfig(; topology="truncated_dilution", matrix_total=1.0))
        r3 = build_reefs(cfg3)
        C3 = build_connectivity(cfg3, r3)
        α, x0 = -log(0.05) / 20_000, 157.57
        kernel(x) = x < x0 ? 0.0 : exp(-α * x)
        trapz(f, a, b; m=200_000) = (h = (b - a) / m; h * (sum(f(a + k * h) for k in 1:(m - 1)) + (f(a) + f(b)) / 2))
        A = r3.area .* r3.k ./ 100
        rad = sqrt.(A ./ π)
        d = 1000 .* SyntheticReefDomain.distance_matrix_km(r3)
        @test all(d[i, i + 1] - rad[i + 1] < x0 for i in 1:2)
        raw = [i == j ? trapz(kernel, 0.0, rad[i]) * A[i] : trapz(kernel, d[i, j] - rad[j], d[i, j]) * A[j]
               for i in 1:3, j in 1:3]
        @test C3 ≈ raw ./ sum(raw) rtol = 1e-4
        @test all(C3[i, i] > 0 for i in 1:3)

        # Links are weighted by the destination: the larger reef receives more
        C2 = td(2, ReefConfig(; n=2, area_m2=[50_000.0, 200_000.0]))
        @test C2[1, 2] > C2[2, 1] > 0

        # One reef: all zeros when it is too small to retain larvae, else the whole total
        @test td(1) == zeros(1, 1)
        @test td(1, ReefConfig(; n=1, area_m2=200_000.0, k_percent=60.0)) ≈ fill(0.00262, 1, 1)

        # Long, weak links fall under ADRIA's 1e-6 cutoff
        C20 = td(20)
        @test 0.0025 < sum(C20) < 0.00262
        @test C20[1, 20] == 0
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
        @test size(D) == (75, 10, 7)  # zero scenario + 6 stochastic
        @test all(D .>= 0)
        @test size(regional[45]) == (75, 7)
        @test mean(D[66:75, :, 2:end]) > mean(D[1:10, :, 2:end])  # warming trend

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
        @test all(Δ[26, 2:end] .≈ 10.0)
        @test all(Δ[setdiff(1:75, 26), :] .== 0.0)
        @test all(reg[45][:, 1] .== 0.0)  # not in the zero scenario

        # Reproducible
        @test build_dhw(cfg, reefs)[1][45] == D
        other = small_cfg(; seed=99)
        @test build_dhw(other, build_reefs(other))[1][45] != D

        # Cyclones: system-wide, group specific, within [0, 1]
        M, intensity = build_cyclones(cfg, 10)
        @test size(M) == (75, 10, 5, 6)  # zero scenario + 5 stochastic
        @test all(0 .<= M .<= 1)
        @test all(M[:, r, :, :] == M[:, 1, :, :] for r in 2:10)
        @test any(intensity .> 0)
        hit = findfirst(>(0), intensity)
        @test issorted(M[hit[1], 1, :, hit[2]]; rev=true)  # vulnerability order

        off = small_cfg(; cyclones=CycloneConfig(; enabled=false))
        @test all(build_cyclones(off, 10)[1] .== 0)

        fixed = small_cfg(; cyclones=CycloneConfig(; annual_prob=0.0, event_years=[2030]))
        _, fi = build_cyclones(fixed, 10)
        @test all(fi[6, 2:end] .> 0)
        @test fi[6, 1] == 0
        @test count(>(0), fi) == 5

        # Waves: ub_med repeated, with as many scenarios as DHW
        W = build_waves(cfg, reefs)
        @test size(W) == size(D)
        @test W[17, :, 3] == reefs.ub_med
    end

    @testset "Zero scenario" begin
        # Scenario 1 is disturbance-free even with frequent, forced and severe events
        cfg = small_cfg(;
            dhw=DHWConfig(;
                baseline=8.0, event_prob=1.0, event_years=[2030], reef_offset=3.0,
                climate=Dict(1 => Dict(), 2 => Dict("trend_per_decade" => 3.0))
            ),
            cyclones=CycloneConfig(; annual_prob=1.0, event_years=[2030])
        )
        synth = build_domain(cfg)
        for label in (1, 2)
            @test size(synth.dhw[label], 3) == 7
            @test all(synth.dhw[label][:, :, ZERO_SCENARIO] .== 0)
            @test all(synth.dhw_regional[label][:, ZERO_SCENARIO] .== 0)
            @test all(synth.dhw[label][:, :, 2:end] .> 0)
        end
        @test size(synth.cyclone_mortality, 4) == 6
        @test all(synth.cyclone_mortality[:, :, :, ZERO_SCENARIO] .== 0)
        @test all(synth.cyclone_intensity[:, ZERO_SCENARIO] .== 0)
        @test all(synth.cyclone_intensity[:, 2:end] .> 0)
        @test size(synth.waves, 3) == 7
        @test synth.waves[:, :, ZERO_SCENARIO] == synth.waves[:, :, 2]

        # ...and in ADRIA, including excluding it from sampling
        path = write_domain(synth, joinpath(TMP, "ZeroScenario"))
        dom = ADRIA.load_domain(path, "2")
        @test all(Array(dom.dhw_scens)[:, :, 1] .== 0)
        @test all(Array(dom.cyclone_mortality_scens)[:, :, :, 1] .== 0)
        @test size(dom.dhw_scens, 3) == 7
        ADRIA.set_factor_bounds!(dom, :dhw_scenario, Tuple(2.0:7.0))
        ADRIA.set_factor_bounds!(dom, :cyclone_mortality_scenario, Tuple(2.0:6.0))
        scens = ADRIA.sample_cf(dom, 32)
        @test all(scens.dhw_scenario .>= 2)
        @test all(scens.cyclone_mortality_scenario .>= 2)
    end

    @testset "Domain README report" begin
        @test SyntheticReefDomain.plots_available()
        cfg = small_cfg(;
            name="Report",
            dhw=DHWConfig(;
                climate=Dict(1 => Dict(), 2 => Dict("trend_per_decade" => 3.0, "event_prob" => 0.3))
            )
        )
        synth = build_domain(cfg)
        path = write_domain(synth, joinpath(TMP, "Report"))
        readme = read(joinpath(path, "README.md"), String)
        figs = ["map.png", "connectivity.png", "cover.png", "dhw_1.png", "dhw_2.png", "cyclones.png"]
        @test sort(readdir(joinpath(path, "figures"))) == sort(figs)
        @test all(filesize(joinpath(path, "figures", f)) > 10_000 for f in figs)
        @test all(occursin("](figures/$f)", readme) for f in figs)
        for heading in ("### Domain", "### Reefs", "### Connectivity", "### Initial cover",
            "### Heat stress (DHW)", "### Cyclones (system-wide)", "### Reefs (resolved values)",
            "## Plots", "## Files")
            @test occursin(heading, readme)
        end
        @test occursin("event_prob = 0.3", readme)  # per-level override listed
        @test count(l -> startswith(l, "| reef_"), split(readme, "
")) == 10
        @test occursin("1 control (no heat stress) + 6 stochastic", readme)

        # Without figures: all settings, no plots, and a note saying how to get them
        path2 = write_domain(synth, joinpath(TMP, "ReportNoFigures"); figures=false)
        readme2 = read(joinpath(path2, "README.md"), String)
        @test !isdir(joinpath(path2, "figures"))
        @test occursin("Plots were not generated", readme2)
        @test occursin("### Reefs (resolved values)", readme2)
        @test !occursin("](figures/", readme2)

        # truncated_dilution lists its own settings instead of self-recruitment and strength
        td_cfg = small_cfg(; n=4, name="ReportTD",
            connectivity=ConnectivityConfig(; topology="truncated_dilution"))
        td_path = write_domain(build_domain(td_cfg), joinpath(TMP, "ReportTD"); figures=false)
        readme3 = read(joinpath(td_path, "README.md"), String)
        @test occursin("| No settlement within | 157.57 m", readme3)
        @test occursin("| Dispersal distance (kernel down to 5%) | 20 km |", readme3)
        @test occursin("| Reefs with local retention | 0 of 4", readme3)
        @test occursin("| Off-diagonal links | 12 |", readme3)
        @test !occursin("| Self-recruitment |", readme3)
        @test occursin("topology=truncated_dilution threshold_m=157.57",
            read(joinpath(td_path, "connectivity", SyntheticReefDomain.connectivity_filename(td_cfg)), String))
    end

    @testset "C~scape export" begin
        cfg = small_cfg(;
            n=3, name="Cs3",
            dhw=DHWConfig(; climate=Dict(45 => Dict(), 85 => Dict("trend_per_decade" => 2.5))),
            cyclones=CycloneConfig(; annual_prob=0.3)
        )
        synth = build_domain(cfg)

        # Options: TOML round trip and validation
        opts = CscapeOptions(; region="Offshore_Northern", connectivity_scale=0.5)
        f = joinpath(TMP, "cscape_export.toml")
        save_cscape_options(f, opts)
        @test load_cscape_options(f) == opts || SyntheticReefDomain.cscape_options_to_dict(load_cscape_options(f)) == SyntheticReefDomain.cscape_options_to_dict(opts)
        @test_throws ArgumentError SyntheticReefDomain.validate(CscapeOptions(; demography_files=["a.rds"]))
        @test_throws ArgumentError SyntheticReefDomain.validate(CscapeOptions(; region=" "))
        @test_throws ArgumentError SyntheticReefDomain.validate(CscapeOptions(; brood_site_value=2.0))

        # Cyclone categories and replicate mapping
        @test cyclone_category.([0.0, 0.2, 0.21, 0.4, 0.6, 0.61, 0.8, 1.0]) == [0, 1, 2, 2, 3, 4, 4, 5]
        @test SyntheticReefDomain.cyclone_scenario_for_replicate.(1:9, 4) == [1, 2, 3, 4, 2, 3, 4, 2, 3]

        inputs = cscape_inputs(synth, opts)
        T, n, R = 75, 3, 7  # replicates = 6 DHW scenarios + control
        @test size(inputs.connectivity) == (n, n, 6)
        @test all(inputs.connectivity[:, :, t] ≈ 0.5 .* permutedims(synth.connectivity) for t in 1:5)
        brood = inputs.connectivity[:, :, 6]
        @test all(brood[i, i] ≈ 0.9 * 0.03 for i in 1:n)
        @test sort(collect(keys(inputs.disturbances))) == [45, 85]
        A = inputs.disturbances[85]
        @test size(A) == (T, n, R, 7)
        @test A[:, :, 2:R, 2] == synth.dhw[85][:, :, 2:R]
        @test all(A[:, :, 1, 2:3] .== 0)               # replicate 1: control
        @test all(in(0:5), A[:, :, :, 3])               # integer categories
        @test any(A[:, :, :, 3] .> 0)
        @test all(A[:, :, :, 4:7] .== 0)                # no COTS
        @test all(A[:, s, :, 3] == A[:, 1, :, 3] for s in 2:n)  # system-wide cyclones
        ic = inputs.init_cover
        @test names(ic) == ["reef_siteid"; CSCAPE_TYPES]
        @test vec(sum(Matrix(ic[:, 2:end]); dims=2)) ≈ vec(sum(synth.cover; dims=1)) .* synth.reefs.k
        @test ic.corym_non_acro == ic.brooder

        # Distance-based connectivity keeps its absolute total in the spawner matrices
        td3 = build_domain(small_cfg(; n=3, name="CsTD",
            connectivity=ConnectivityConfig(; topology="truncated_dilution")))
        @test sum(cscape_inputs(td3, CscapeOptions()).connectivity[:, :, 1]) ≈ 0.00262

        # Writing without R: CSV inputs, init cover, script; force keeps added files
        dir = joinpath(TMP, "Cs3_cscape")
        r = write_cscape_project(synth, dir; options=opts, run_r=false)
        @test !r.built
        for fname in ("build_cscape_project.R", "README.md", "cscape_export.toml",
            "generator_config.toml", joinpath("data", "initialisationData.csv"),
            joinpath("inputs", "sites.csv"), joinpath("inputs", "connectivity.csv"),
            joinpath("inputs", "scenarios.csv"), joinpath("inputs", "disturbances_45.csv"),
            joinpath("inputs", "disturbances_85.csv"))
            @test isfile(joinpath(dir, fname))
        end
        @test countlines(joinpath(dir, "inputs", "disturbances_45.csv")) == 1 + T * n * R
        @test countlines(joinpath(dir, "inputs", "scenarios.csv")) == 3
        @test occursin("demog_inputs_brooder.rds", read(joinpath(dir, "README.md"), String))
        @test_throws ArgumentError write_cscape_project(synth, dir; options=opts, run_r=false)
        touch(joinpath(dir, "data", "demog_inputs_acro_table.rds"))
        write_cscape_project(synth, dir; options=opts, run_r=false, force=true)
        @test isfile(joinpath(dir, "data", "demog_inputs_acro_table.rds"))
        other = mkpath(joinpath(TMP, "NotAnExport"))
        touch(joinpath(other, "keep.txt"))
        @test_throws ArgumentError write_cscape_project(synth, other; run_r=false, force=true)

        # Zip downloads
        Z = SyntheticReefDomain.ZipArchives
        zd = Z.ZipReader(domain_zip(synth; figures=false))
        @test "Cs3/datapackage.json" in Z.zip_names(zd)
        @test "Cs3/DHWs/dhwRCP85.nc" in Z.zip_names(zd)
        zc = Z.ZipReader(cscape_zip(synth, opts; run_r=false))
        @test "Cs3_cscape/build_cscape_project.R" in Z.zip_names(zc)
        @test "Cs3_cscape/inputs/disturbances_85.csv" in Z.zip_names(zc)

        # With R and cscape installed: build with R, then let C~scape read it (dry run) and
        # simulate one replicate using the dummy demography of CScape.jl's 00_Verify project
        status = rscript_status()
        demog = joinpath(@__DIR__, "..", "..", "CScape.jl", "sandbox", "projects", "00_Verify", "data")
        if status.rscript !== nothing && status.cscape && isdir(demog)
            rdir = joinpath(TMP, "Cs3_cscape_R")
            rr = write_cscape_project(synth, rdir; options=CscapeOptions())
            @test rr.built
            for fname in ("ScenarioID.xlsx", joinpath("data", "reef_sites_Cs3.rds"),
                joinpath("data", "connectivity_Cs3.rds"), joinpath("data", "temporal_Cs3_45.rds"),
                joinpath("data", "temporal_Cs3_85.rds"))
                @test isfile(joinpath(rdir, fname))
            end
            for fname in readdir(demog)
                startswith(fname, "demog_inputs_") && cp(joinpath(demog, fname), joinpath(rdir, "data", fname))
            end
            check = joinpath(TMP, "check_cscape.R")
            write(check, """
            fpath <- commandArgs(trailingOnly = TRUE)[1]
            suppressMessages(library(cscape))
            tmp <- readRDS(file.path(fpath, "data", "temporal_Cs3_85.rds"))
            stopifnot(identical(dim(tmp\$disturbances), c(75L, 3L, 7L, 7L)))
            stopifnot(all(tmp\$disturbances[, , 1, 2:3] == 0))
            invisible(run_cscape(85, fpath, rep_index = 2, dry_run = TRUE))
            invisible(run_cscape(85, fpath, rep_index = 2))
            stopifnot(file.exists(file.path(fpath, "model_outputs", "Array_scenario_85_rep_2_draw_NA.rds")))
            cat("CSCAPE_OK\n")
            """)
            out = read(pipeline(`$(status.rscript) $check $rdir`; stderr=devnull), String)
            @test occursin("CSCAPE_OK", out)
        else
            @info "Skipping the C~scape R round trip (needs Rscript, the cscape package and CScape.jl's 00_Verify demography)."
        end
    end

    @testset "ADRIA round trip" begin
        cases = [(n, "distance_decay") for n in (1, 2, 20)]
        append!(cases, [(10, t) for t in TOPOLOGIES])
        for (n, topology) in cases
            name = "RT_$(n)_$(topology)"
            cfg = small_cfg(; n=n, name=name, connectivity=ConnectivityConfig(; topology=topology))
            synth = build_domain(cfg)
            path = write_domain(synth, joinpath(TMP, name); figures=false)
            @test isfile(joinpath(path, "spatial", "RT_$(n)_$(topology).gpkg"))

            dom = ADRIA.load_domain(path, "45")
            ids = synth.reefs.reef_siteid
            @test dom.loc_ids == ids
            @test dom.loc_data.reef_siteid == ids
            @test Array(dom.conn) == synth.connectivity
            @test Array(dom.dhw_scens) == synth.dhw[45]
            @test Array(dom.wave_scens) == synth.waves
            @test size(dom.cyclone_mortality_scens) == (75, n, 5, 6)
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
            @test write_domain(synth, path; force=true, figures=false) == path
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
            for (n, topology) in ((2, "distance_decay"), (10, "distance_decay"),
                (20, "distance_decay"), (10, "truncated_dilution"))
                name = "Run_$(n)_$(topology)"
                cfg = small_cfg(; n=n, name=name, connectivity=ConnectivityConfig(; topology=topology))
                path = write_domain(build_domain(cfg), joinpath(TMP, name); figures=false)
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
