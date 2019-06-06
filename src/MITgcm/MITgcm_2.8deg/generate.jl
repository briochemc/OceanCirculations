#=
This code is to load the MITgcm matrix and grid
and adapt it to AIBECS.

Importantly, the TMM uses a different indexing of the 3D grid.
In TMM, the lons are traversed first
=#

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using MAT                   # For loading MAT format
using BSON                  # For saving circulation as BSON format
using Unitful, UnitfulAstro # for units
using OceanGrids            # To store the grid

# No download, this is a local test file

println("Reading grid file (grid data)")
mat_boxes_file = "/Users/benoitpasquier/Data/MITgcm_2.8deg/Matrix5/Data/boxes.mat"
boxes = matread(mat_boxes_file)
mat_grd_file = "/Users/benoitpasquier/Data/MITgcm_2.8deg/grid.mat"
grd = matread(mat_grd_file)

println("Reading Matrix file (circulation matrix)")
mat_Matrix_file = "/Users/benoitpasquier/Data/MITgcm_2.8deg/Matrix5/TMs/matrix_nocorrection_annualmean.mat"
vars = matread(mat_Matrix_file)

println("  Circulation (T)")
T = -vars["Aexpms"] * upreferred(u"1/yr") # TODO check unit and direction is correct!

println("  Grid")
lat = vec(grd["y"]) * u"°"
lon = vec(grd["x"]) * u"°"
depth = vec(grd["z"]) * u"m"
nlat, nlon, ndepth = Int(grd["ny"]), Int(grd["nx"]), Int(grd["nz"])
δx_3D = repeat(reshape(permutedims(grd["dx"], [2,1]), (nlat,nlon,1)), outer=(1,1,ndepth)) * u"m"
δy_3D = repeat(reshape(permutedims(grd["dy"], [2,1]), (nlat,nlon,1)), outer=(1,1,ndepth)) * u"m"
δz_3D = repeat(reshape(grd["dznom"], (1,1,ndepth)), outer=(nlat,nlon,1)) * u"m"
volume_3D = δx_3D .* δy_3D .* δz_3D

wet3d_orig = convert(BitArray{3}, grd["bathy"])
iwet_orig = findall(vec(wet3d_orig))
wet3d = permutedims(wet3d_orig, [2, 1, 3])
iwet = findall(vec(wet3d))

permuted_iwet = permutedims(LinearIndices(size(wet3d)), [2, 1, 3])[iwet_orig]
p = sortperm(permuted_iwet)

# TODO here figure out the index conversions

volume_3D = δx_3D .* δy_3D .* δz_3D
depth_top_3D = grd["ZU"] * u"m"
depth_top = depth_top_3D[1,1,:]
A_2D = grd["TAREA"] * u"m^2"
nlat, nlon, ndepth = size(grd["TLONG"])
nboxes = length(grd["TLONG"])
OceanCurvilinearGrid(
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

bson_dir = joinpath(splitpath(mat_MTM0_file)[1:end-2]... ,"POPX3_BSON")
bson_file = joinpath(bson_dir, "POPX3.bson")
isdir(bson_dir) || mkdir(bson_dir)
isfile(bson_file) && rm(bson_file)
BSON.@save bson_file grid wet3D T
