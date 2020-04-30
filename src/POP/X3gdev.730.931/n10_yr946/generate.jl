#=
The code below loads the grid information and transport matrix
from the POP ocean circulation model (from Ann Bardin)
and adapt it to AIBECS.
=#

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using MAT                   # For loading MAT format
using BSON                  # For saving circulation as BSON format
using Unitful               # for units
using OceanGrids            # To store the grid

# No download, this is a local test file

println("Reading MET file (grid data)")
data_path = "/Users/benoitpasquier/Data"
mat_MET_file = joinpath(data_path, "abardin/POPOFFLINE_ops/X3gdev.730.931/n10_yr946/MET.mat")
grd = matread(mat_MET_file)["MET"]

println("Reading MTM0 file (circulation matrix)")
mat_MTM0_file = joinpath(data_path, "abardin/POPOFFLINE_ops/X3gdev.730.931/n10_yr946/MTM0.mat")
vars = matread(mat_MTM0_file)

println("  Circulation (T = A + H + D)")
T = -vars["T"] * u"1/s" # TODO check unit and direction is correct!
A = -vars["A"] * u"1/s" # TODO check unit and direction is correct!
H = -vars["H"] * u"1/s" # TODO check unit and direction is correct!
D = -vars["D"] * u"1/s" # TODO check unit and direction is correct!

println("  Grid")
lat_3D = grd["TLONG"] * u"°"
lon_3D = grd["TLAT"] * u"°"
depth_3D = grd["ZT"] * u"m"
δx_3D = grd["DXT"] * u"m"
δy_3D = grd["DYT"] * u"m"
δz_3D = grd["DZT"] * u"m"
volume_3D = δx_3D .* δy_3D .* δz_3D
depth_top_3D = grd["ZU"] * u"m"
depth_top = depth_top_3D[1,1,:]
A_2D = grd["TAREA"] * u"m^2"
nlat, nlon, ndepth = size(grd["TLONG"])
nboxes = length(grd["TLONG"])
grid = OceanCurvilinearGrid(
                            lat_3D,
                            lon_3D,
                            depth_3D,
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
wet3D = falses(nlat, nlon, ndepth)
wet3D[Int.(vec(grd["iocn"]))] .= true

bson_dir = joinpath(data_path, "OceanGrids")
println("Saving as BSON file in $bson_dir")
bson_file = joinpath(bson_dir, "POPX3.bson")
isdir(bson_dir) || mkdir(bson_dir)
isfile(bson_file) && rm(bson_file)
BSON.@save bson_file grid wet3D T
