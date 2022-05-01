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

# Create registry entry for OCIM1 MAT file
println("Registering OCIM1 MAT data (making it available locally)")
register(
    DataDep(
        "OCIM1_MATLAB",
        """
        References:
        - DeVries, T. (2014), The oceanic anthropogenic CO2 sink: Storage, air‐sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.
        - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
        """,
        "https://dl.dropboxusercontent.com/s/n0j1swt4bqgbj75/CTL.mat",
        "347a6bb17676f957c6f56c46313d11298e4764f6f5d0304b7ea8374108fb442f",
        fetch_method = fallback_download
    )
)

println("Reading OCIM1 MAT file")
mat_file = @datadep_str string("OCIM1_MATLAB/", "CTL.mat")
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

data_path = "/Users/benoitpasquier/Data"
jld2_dir = joinpath(data_path, "OceanGrids")
println("Saving as JLD2 file in $jld2_dir")
jld2_file = joinpath(jld2_dir, "OCIM1_CTL.jld2")
isdir(jld2_dir) || mkdir(jld2_dir)
isfile(jld2_file) && rm(jld2_file)
jldsave(jld2_file, true; grid, T)
