#=
The code below loads the grid information and transport matrix
from the OCCA output by Gael Forget
and adapts them to AIBECS.
=#

using SparseArrays          # For sparse matrix in OCIM
using MAT                   # For loading OCIM in MAT format
using BSON                  # For saving circulation as BSON format
using Unitful, UnitfulAstro # for units
using OceanGrids            # To store the grid

println("Reading OCCA MAT file")
data_path = "/Users/benoitpasquier/Data"
mat_file = joinpath(data_path, "OCCA-TRM", "matrix_10.mat")
vars = matread(mat_file)


println("  Wet boxes")
wet3D_orig = vars["mask2x2"] .== 1
iwet_orig = findall(vec(wet3D_orig))
wet3D = permutedims(wet3D_orig, [2,1,3]) # swap lat and lon
iwet = findall(vec(wet3D))
permuted_iwet = permutedims(LinearIndices(size(wet3D)), [2, 1, 3])[iwet_orig]
p = sortperm(permuted_iwet)

println("  Circulation")
T = (-vars["A_adv"] - vars["A_dif"])[p,p] * u"1/s" # `p` swaps lat and lon? TODO check

println("  Grid")
lat = vars["YC2x2"][1,:] * u"°"
lon = vars["XC2x2"][:,1] * u"°"
depth = vec(vars["RC2x2"]) * u"m"
δlat = fill(2u"°", size(lat))
δlon = fill(2u"°", size(lon))
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
