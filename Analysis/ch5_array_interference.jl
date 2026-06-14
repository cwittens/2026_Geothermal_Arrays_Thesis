using Pkg
Pkg.activate(@__DIR__)

using JLD2: load
using Statistics: mean
using Plots
using LaTeXStrings
using GeothermalWells

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()

data_dir_2x2 = joinpath(simulation_dir(), "array_studies_2x2")
data_dir_3x3 = joinpath(simulation_dir(), "array_studies_3x3")


default(
    grid=true,
    box=:on,
    size=(600, 450),
    dpi=300,
    titlefont=font(16),
    linewidth=3, gridlinewidth=2,
    markersize=4, markerstrokewidth=0.1,
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=16, yguidefontsize=16,
    ztickfontsize=14, zguidefontsize=16,
    legendfontsize=14,
    fontfamily="Computer modern"
)


spacings = collect(0:10:120)     # metres; idx 1 = 0 (single well), idx i = (i-1)*10
nfiles   = length(spacings)

α = 2.55 / 2.356e6               # rock thermal diffusivity [m^2/s]
yr = 3600 * 24 * 365

# outlet (z≈0) and well bottom (z≈h) temperature time series for one borehole
function well_temps(sv, cc, bh)
    ij = [(i, j) for (i, j, b) in cc.Idx_list_Outer if b == bh]
    gridz = cc.gridz
    h = cc.boreholes[bh].h
    k_top = 1
    k_bot = findlast(gridz .<= h)
    T_top = [mean(T[k_top, j, i] for (i, j) in ij) for T in sv.saveval]
    T_bot = [mean(T[k_bot, j, i] for (i, j) in ij) for T in sv.saveval]
    return T_top, T_bot
end

# neighbour count from coordinates: center=4, edge=3, corner=2
function well_class(bh, cc; tol=1e-6)
    xc, yc = cc.boreholes[bh].xc, cc.boreholes[bh].yc
    n = (abs(xc) > tol) + (abs(yc) > tol)
    return n == 0 ? :center : (n == 1 ? :edge : :corner)
end


function extract_2x2(path)
    sv = load(path, "saved_values")
    cc = load(path, "cache_cpu")
    times = sv.t ./ yr
    T_top, T_bot = well_temps(sv, cc, 1)   # all 4 wells equivalent by symmetry
    return (; times, T_top, T_bot)
end

function extract_3x3(path)
    sv = load(path, "saved_values")
    cc = load(path, "cache_cpu")
    times = sv.t ./ yr
    reps = Dict{Symbol,Int}()
    for bh in eachindex(cc.boreholes)
        get!(reps, well_class(bh, cc), bh)   # one representative per class
    end
    classes = Dict(cls => well_temps(sv, cc, bh) for (cls, bh) in reps)
    return (; times, classes)
end

data_2x2 = Vector{Any}(undef, nfiles)
for i in 1:nfiles
    path = joinpath(data_dir_2x2, "array_study_$(lpad(i,2,'0')).jld2")
    isfile(path) || (@warn "missing 2x2 idx $i"; continue)
    data_2x2[i] = extract_2x2(path)
    GC.gc()
    @info "2x2 done: spacing $(spacings[i]) m"
end

data_3x3 = Vector{Any}(undef, nfiles)
for i in 1:nfiles
    path = i == 1 ? joinpath(data_dir_2x2, "array_study_01.jld2") :
                    joinpath(data_dir_3x3, "3x3_array_study_$(lpad(i,2,'0')).jld2")
    isfile(path) || (@warn "missing 3x3 idx $i"; continue)
    data_3x3[i] = extract_3x3(path)
    GC.gc()
    @info "3x3 done: spacing $(spacings[i]) m"
end

# =============================================================================
# Part A
# =============================================================================
begin
T_out_final  = [data_2x2[i].T_top[end] for i in 1:nfiles]
T_turn_final = [data_2x2[i].T_bot[end] for i in 1:nfiles]

order = [2:nfiles; 1]            # spacing 0 (single well) plotted as the ∞ limit
x = 1:nfiles

pA = plot(x, T_out_final[order], label="Outlet (z=0)", marker=:circle,
          legend=:bottomright,
           ylims=(19.5, 50.5),
           xtickfontsize=13, ytickfontsize=13,
           )
plot!(pA, x, T_turn_final[order], label="Well bottom (z=3000 m)", marker=:square)
xticks!(pA, x, [string.(10:10:120); L"\infty"])
xlabel!(pA, "Borehole spacing [m]")
ylabel!(pA, "Temperature [°C]")
# title!(pA, "2x2 array, 20 years")

# todo make // lines on x axis to indicate break between finite spacings and single-well limit

@info savefig(pA, joinpath(plots_dir(), "array_2x2_static.pdf"))
end
# =============================================================================
# Part B
# =============================================================================
begin
times  = data_2x2[1].times
ntimes = length(times)

# outlet temperature matrix: [spacing, time]
Tmat = fill(NaN, nfiles, ntimes)
for i in 1:nfiles
    Tmat[i, :] = data_2x2[i].T_top
end
T_single = Tmat[1, :]                       # spacing 0 reference

pos = 2:nfiles                              # spacings 10..120 m
pos_s = spacings[pos]


# B: characteristic interference distance s_c(t), test s_c ∝ √t
function crossing(s, dev, thr)              # dev decreasing in s; first crossing
    for k in 1:length(s)-1
        if (dev[k] - thr) * (dev[k+1] - thr) <= 0 # one below one above threshold crossing
            f = (thr - dev[k]) / (dev[k+1] - dev[k])
            return s[k] + f * (s[k+1] - s[k]) # interpolate between s[k] and s[k+1]
        end
    end
    return NaN
end

thr = 0.5
s_c, t_c = Float64[], Float64[]
for tk in 1:ntimes
    dev = T_single[tk] .- Tmat[pos, tk]
    sc = crossing(pos_s, dev, thr)
    isnan(sc) || (push!(s_c, sc); push!(t_c, times[tk]))
end

γ = sum(sqrt.(t_c) .* s_c) / sum(t_c)       # LS fit of s_c = γ√t (through origin)
t_dense = range(minimum(t_c), maximum(t_c), 100)
γ / sqrt(2.55 / 2.356e6 * 365*24*3600)

pB = scatter(t_c, s_c, 
label = "Onset of interference",
xticks = [1, 5, 10, 15, 20],
marker=:o, markersize=5, color=1, legend=:topleft
)

plot!(pB, t_dense, γ .* sqrt.(t_dense), label="∝ √t fit")
xlabel!(pB, "Time [yr]"); ylabel!(pB,"Interference distance " * L"\textrm{s}_\textrm{c}" * " [m]")


@info savefig(pB, joinpath(plots_dir(),  "interference_sqrt_t.pdf"))
end



# =============================================================================
# Part C
# =============================================================================

begin
final_T(d, cls) = haskey(d.classes, cls) ? d.classes[cls][1][end] : NaN  # [1]=T_top

T_corner = [isassigned(data_3x3, i) ? final_T(data_3x3[i], :corner) : NaN for i in 1:nfiles]
T_edge   = [isassigned(data_3x3, i) ? final_T(data_3x3[i], :edge)   : NaN for i in 1:nfiles]
T_center = [isassigned(data_3x3, i) ? final_T(data_3x3[i], :center) : NaN for i in 1:nfiles]

single = data_3x3[1].classes[:center][1][end]   # single-well reference (= ∞ limit)

pC = plot(legend=:bottomright, ylims=(0, 36))
plot!(pC, x, [T_corner[2:nfiles]; single], marker=:circle,  label="Corner Well")
plot!(pC, x, [T_edge[2:nfiles];   single], marker=:square,  label="Edge Well")
plot!(pC, x, [T_center[2:nfiles]; single], marker=:diamond, label="Central Well")
xticks!(pC, x, [string.(10:10:120); L"\infty"])
xlabel!(pC, "Borehole spacing [m]"); ylabel!(pC, "Outlet temperature [°C]")
# title!(pC, "3×3 array by well position, 20 years")

@info savefig(pC, joinpath(plots_dir(), "array_3x3_classes.pdf"))

end