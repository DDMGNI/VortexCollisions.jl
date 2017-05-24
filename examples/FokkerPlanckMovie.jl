
#runid = "FokkerPlanck_dt1E-6_nt10000"
#runid = "FokkerPlanck_exp_cos_dt1E-6_nt10000"
#runid = "FokkerPlanck_sin4xsin2y_dt1E-6_nt10000"
#runid = "FokkerPlanck_sinx4_dt1E+3_nt100"
#runid = "FokkerPlanck_sinx4siny4_dt1E-3_nt10000"
runid = "FokkerPlanck_sinx4sinyetc_dt1E-3_nt10000"


using PyCall

@pyimport matplotlib as mpl
mpl.use("Agg")
mpl.interactive(false)

@pyimport matplotlib.gridspec as gspec
@pyimport matplotlib.colors   as clrs

using PyPlot
using HDF5


h5 = h5open(runid * ".h5", "r")

nplot = size(h5["ω"],3)
#nplot = 200
#nplot = 2

xaxis = collect(1:nplot)-1;

ϕ = view(h5["ϕ"][:,:,1], :, :, 1)
ω = view(h5["ω"][:,:,1], :, :, 1)

ϕ_min = minimum(ϕ)
ϕ_max = maximum(ϕ)
ϕ_dif = ϕ_max - ϕ_min

ω_min = minimum(ω)
ω_max = maximum(ω)
ω_dif = ω_max - ω_min


h = zeros(nplot)
for n in 1:nplot
    h[n] = sum_kbn(h5["ω"][:,:,n] .* h5["ϕ"][:,:,n])
end

h_err = (h-h[1])/h[1]
h_max = maximum(h_err)
h_min = minimum(h_err)

ℰ = zeros(nplot)
for n in 1:nplot
    ℰ[n] = sum_kbn(h5["ω"][:,:,n] .* h5["ω"][:,:,n])
end

ℰ_err = (ℰ-ℰ[1])/ℰ[1]
ℰ_max = maximum(ℰ_err)
ℰ_min = minimum(ℰ_err)


fig = figure(figsize=(16,10))

subplots_adjust(left=0.05, right=0.98, top=0.96, bottom=0.07)

gs = gspec.GridSpec(2,1, height_ratios=(2,1))

gs0 = gspec.GridSpecFromSubplotSpec(1,3, subplot_spec=get(gs, 0,0))
gs1 = gspec.GridSpecFromSubplotSpec(1,2, subplot_spec=get(gs, 1,0))

axϕ = subplot(get(gs0, 0,0))
axω = subplot(get(gs0, 1,0))
axs = subplot(get(gs0, 2,0))
axh = subplot(get(gs1, 0,0))
axℰ = subplot(get(gs1, 1,0))


ncnt = 40

#ϕnorm = clrs.Normalize(vmin=ϕ_min - 0.1 * ϕ_dif, vmax=ϕ_max + 0.1 * ϕ_dif)
#ωnorm = clrs.Normalize(vmin=ω_min - 0.1 * ω_dif, vmax=ω_max + 0.1 * ω_dif)

ϕlevels = linspace(ϕ_min - ϕ_dif, ϕ_max + ϕ_dif, ncnt)
ωlevels = linspace(ω_min - ω_dif, ω_max + ω_dif, ncnt)

pcmϕ = axϕ[:pcolormesh](ϕ)
pcmω = axω[:pcolormesh](ω)

#cntϕ = axϕ[:contour](ϕ, ncnt, colors="k", norm=ϕnorm)
#cntω = axω[:contour](ω, ncnt, colors="k", norm=ωnorm)
cntϕ = axϕ[:contour](ϕ, ϕlevels, colors="k")
cntω = axω[:contour](ω, ωlevels, colors="k")

scat = axs[:scatter](ϕ[:], ω[:])

linesh = axh[:plot](xaxis[1:1], h_err[1:1])[1]
linesℰ = axℰ[:plot](xaxis[1:1], ℰ_err[1:1])[1]

axϕ[:set_title]("ϕ",         fontsize=18)
axω[:set_title]("ω",         fontsize=18)
axs[:set_title]("ω(ϕ)",      fontsize=18)
axh[:set_title]("Energy",    fontsize=18)
axℰ[:set_title]("Enstrophy", fontsize=18)

#axϕ[:colorbar]()
#axω[:colorbar]()

axs[:set_xlabel]("ϕ",   fontsize=16)
axs[:set_ylabel]("ω",   fontsize=16)

axh[:set_xlim]((xaxis[1], xaxis[end]))
axh[:set_ylim]((h_min, h_max))
axh[:set_xlabel]("t",         fontsize=16)
axh[:set_ylabel]("(h-h₀)/h₀", fontsize=16)

axℰ[:set_xlim]((xaxis[1], xaxis[end]))
axℰ[:set_ylim]((ℰ_min, ℰ_max))
axℰ[:set_xlabel]("t",         fontsize=16)
axℰ[:set_ylabel]("(ℰ-ℰ₀)/ℰ₀", fontsize=16)


for n in 1:nplot
	
	println("Plotting frame " * @sprintf("%8i", n-1))
    
	ω = view(h5["ω"][:,:,n], :, :, 1)
	ϕ = view(h5["ϕ"][:,:,n], :, :, 1)

    for coll in cntϕ[:collections]
		coll[:remove]()
    end
        
	for coll in cntω[:collections]
		coll[:remove]()
	end
    
    pcmϕ[:set_array](ϕ[:])
    pcmω[:set_array](ω[:])
	
#	cntϕ = axϕ[:contour](ϕ, ncnt, colors="k", norm=ϕnorm)
#	cntω = axω[:contour](ω, ncnt, colors="k", norm=ωnorm)
	cntϕ = axϕ[:contour](ϕ, ϕlevels, colors="k")
	cntω = axω[:contour](ω, ωlevels, colors="k")

    scat[:set_offsets](hcat(ϕ[:], ω[:]))
	
    linesh[:set_data](xaxis[1:n], h_err[1:n])
    linesℰ[:set_data](xaxis[1:n], ℰ_err[1:n])

	savefig(runid * "_movie_" * @sprintf("%08i", n-1) * ".png", dpi=150)

end


close(fig)

close(h5)
