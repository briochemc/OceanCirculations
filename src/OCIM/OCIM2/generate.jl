#=
The code below loads the grid information and transport matrix
from the OCIM (from Tim DeVries)
and adapt it to AIBECS.
=#

using SparseArrays          # For sparse matrix in OCIM
using DataDeps              # For storage location of data
using Downloads
ENV["DATADEPS_ALWAYS_ACCEPT"] = true # so it does not ask me
using MAT                   # For loading OCIM in MAT format
using JLD2                  # For saving circulation as JLD2 format
using Unitful               # for units
using OceanGrids            # To store the grid

# fallback for download
function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end

url_start = "dl.dropboxusercontent.com/sh/eun1dkwjeqvax80/"

OCIM2_versions = Dict(
    :OCIM2_CTL_He             => "$url_start/AADDlMLcnio1Fnz3ifL8ih97a/OCIM2_CTL_He.mat",
    :OCIM2_CTL_noHe           => "$url_start/AADSHmqj7kLoXAqEy5iwNjQra/OCIM2_CTL_noHe.mat",
    :OCIM2_KiHIGH_He          => "$url_start/AAA6IM3Z_hWrjs1JjyBTZU1da/OCIM2_KiHIGH_He.mat",
    :OCIM2_KiHIGH_noHe        => "$url_start/AAAwAR-v4PJb0gi9-YRdJjABa/OCIM2_KiHIGH_noHe.mat",
    :OCIM2_KiLOW_He           => "$url_start/AAD4Q5Ney1SP6DBpvBtMx31La/OCIM2_KiLOW_He.mat",
    :OCIM2_KiLOW_noHe         => "$url_start/AAAspy9g9nguTJH6W1WeGmPGa/OCIM2_KiLOW_noHe.mat",
    :OCIM2_KvHIGH_He          => "$url_start/AADI2g7bze8GTC0LV6r3ovHra/OCIM2_KvHIGH_He.mat",
    :OCIM2_KvHIGH_KiHIGH_noHe => "$url_start/AAB6cH080CMDEBNefTGZaVJsa/OCIM2_KvHIGH_KiHIGH_noHe.mat",
    :OCIM2_KvHIGH_KiLOW_He    => "$url_start/AAAwNU9F_hkPv4DFTAgQM1hKa/OCIM2_KvHIGH_KiLOW_He.mat",
    :OCIM2_KvHIGH_KiLOW_noHe  => "$url_start/AABxfubJtPy2JdqMTRF95c5Ha/OCIM2_KvHIGH_KiLOW_noHe.mat",
    :OCIM2_KvHIGH_noHe        => "$url_start/AADPh-3Ix3-pThNwYBM4x5XNa/OCIM2_KvHIGH_noHe.mat"
)


for (name,url) in OCIM2_versions
    # Create registry entry for OCIM2 MAT file
    println("Registering $name")
    register(
        DataDep(
            "$(name)_MATLAB",
            """
            References:
            - DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716
            """,
            url,
            sha2_256,
            fetch_method = fallback_download
        )
    )

    println("Reading $name.mat file")
    mat_file = @datadep_str string("$(name)_MATLAB/$name.mat")
    vars = matread(mat_file)


    println("  Circulation")
    T = -vars["output"]["TR"] * upreferred(u"1/yr")

    println("  Wet boxes")
    wet3D = convert(BitArray{3}, vars["output"]["M3d"])

    println("  Grid")
    grd = vars["output"]["grid"]
    lat = vec(grd["yt"]) * u"°"
    lon = vec(grd["xt"]) * u"°"
    depth = vec(grd["zt"]) * u"m"
    δlat = vec(grd["dyt"]) * u"°"
    δlon = vec(grd["dxt"]) * u"°"
    δdepth = vec(grd["dzt"]) * u"m"
    lat_3D = grd["YT3d"] * u"°"
    lon_3D = grd["XT3d"] * u"°"
    depth_3D = grd["ZT3d"] * u"m"
    δy = grd["DYT3d"][:,1,1] * u"m"
    δx_3D = grd["DXT3d"] * u"m"
    δy_3D = grd["DYT3d"] * u"m"
    δz_3D = grd["DZT3d"] * u"m"
    volume_3D = δx_3D .* δy_3D .* δz_3D
    depth_top_3D = grd["ZW3d"] * u"m"
    depth_top = depth_top_3D[1,1,:]
    A_2D = grd["Areat"] * u"m^2"
    nlat, nlon, ndepth = size(grd["XT3d"])
    nboxes = length(grd["XT3d"])
    grid = OceanRectilinearGrid(
                     lat,
                     lon,
                     depth,
                     δlat,
                     δlon,
                     δdepth,
                     lat_3D,
                     lon_3D,
                     depth_3D,
                     δy,
                     δx_3D,
                     δy_3D,
                     δz_3D,
                     volume_3D,
                     depth_top,
                     depth_top_3D,
                     A_2D,
                     wet3D,
                     nlon,
                     nlat,
                     ndepth,
                     nboxes
                    )

    println("  He fluxes")
    Base.vec(nothing) = nothing
    Base.:*(::Nothing, ::Unitful.FreeUnits) = nothing
    He3Flux = vec(get(vars["output"], "J_He3", nothing)) * u"mol/m^3/yr" # vector of ³He flux
    He4Flux = vec(get(vars["output"], "J_He4", nothing)) * u"mol/m^3/yr" # vector of ⁴He flux

    data_path = "/Users/benoitpasquier/Data"
    jld2_dir = joinpath(data_path, "OceanGrids")
    println("Saving as JLD2 file in $jld2_dir")
    jld2_file = joinpath(jld2_dir, "$name.jld2")
    isdir(jld2_dir) || mkdir(jld2_dir)
    isfile(jld2_file) && rm(jld2_file)

    jldsave(jld2_file, true; grid, T, He3Flux, He4Flux)

end

