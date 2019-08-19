#=
The code below loads the grid information and transport matrix
from the OCIM0 (from Francois Primeau and Tim DeVries)
and adapt them to AIBECS.
=#

using SparseArrays          # For sparse matrix in OCIM
using MAT                   # For loading OCIM in MAT format
using BSON                  # For saving circulation as BSON format
using Unitful, UnitfulAstro # for units
using OceanGrids            # To store the grid



println("Reading OCIM0 MAT file")
data_path = "/Users/benoitpasquier/Data"
mat_file = joinpath(data_path, "OCIM0.1.mat")
vars = matread(mat_file)

println("  Circulation")
T = vars["T"] * u"1/s"

println("  Grid")
grd = vars["grid"]
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
                 nlon,
                 nlat,
                 ndepth,
                 nboxes
                )

println("  Wet boxes")
wet3D = convert(BitArray{3}, vars["M3d"])

data_path = "/Users/benoitpasquier/Data"
bson_dir = joinpath(data_path, "OceanGrids")
println("Saving as BSON file in $bson_dir")
bson_file = joinpath(bson_dir, "OCIM0.1.bson")
isdir(bson_dir) || mkdir(bson_dir)
isfile(bson_file) && rm(bson_file)
BSON.@save bson_file grid wet3D T
