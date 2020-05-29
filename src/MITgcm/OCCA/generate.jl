#=
The code below loads the grid information and transport matrix
from the OCCA output by Gael Forget
and adapts them to AIBECS.
=#

using SparseArrays          # For sparse matrix
using LinearAlgebra
using MAT                   # For loading OCCA in MAT format
using BSON                  # For saving circulation as BSON format
using Unitful               # for units
using Unitful: °, m, km, s  # for units
using OceanGrids            # To store the grid
using Distances             # to compute distances from (lat,lon)

println("Reading OCCA MAT file")
data_path = "/Users/benoitpasquier/Data"
# Note: The OCCA matrix by Gael Forget was sent to me in an email
# So it is not available in a public repository (yet?)
mat_file = joinpath(data_path, "OCCA-TRM", "matrix_10.mat")
vars = matread(mat_file)


println("  Wet boxes")
wet3D_orig = vars["mask2x2"] .== 1
wet3D = permutedims(wet3D_orig, (2,1,3)) # swap lat and lon
iwet = findall(vec(wet3D))

"""
    edges(e1, x)

Compute the edges given the center points `x` and first edge `e1`.

Basically assume that the `x` values are centered, so edges `e` must be at 
equal distances from the centers. I.e., starting from a `x` that looks like
  x---x----x-------x
create the `e` such that they look like
e-x-e-x-e--x--e----x----e
"""
function edges(e1, x::Vector{T}) where T
    e = Vector{T}(undef, length(x)+1)
    e[1] = e1
    for i in eachindex(x)
        e[i+1] = 2x[i] - e[i]
    end
    return e
end

println("  Grid")
lat = vars["YC2x2"][1,:]°
lon = vars["XC2x2"][:,1]°
depth = vec(vars["RC2x2"])m
elat = edges(-80.0°, lat)
elon = edges(0.0°, lon)
edepth = edges(0m, depth)
δlat = diff(elat)
δlon = diff(elon)
δdepth = diff(edepth)
nlat, nlon, ndepth = length.((lat, lon, depth))
lat_3D = repeat(lat, outer=(1, nlon, ndepth))
lon_3D = repeat(reshape(lon, (1, nlon, 1)), outer=(nlat, 1, ndepth))
depth_3D = repeat(reshape(depth, (1, 1, ndepth)), outer=(nlat, nlon, 1))
earthdistance(p1,p2) = Haversine(6371e3)(ustrip.(p1), ustrip.(p2))m
δy = [earthdistance((0.0°, elat[i+1]), (0.0°, elat[i])) for i in eachindex(lat)]
δy_2D = [earthdistance((lon, elat[i+1]), (lon, elat[i])) for i in eachindex(lat), lon in lon]
δx_2D = [earthdistance((elon[i+1], lat), (elon[i], lat)) for lat in lat, i in eachindex(lon)]
δx_3D = repeat(δx_2D, outer=(1,1,ndepth))
δy_3D = repeat(δy_2D, outer=(1,1,ndepth))
δz_3D = repeat(reshape(δdepth, (1,1,ndepth)), outer=(nlat,nlon,1))
volume_3D = permutedims(vars["vol2x2"], (2,1,3))m^3 # swap lat and lon
depth_top = edepth[1:end-1]
depth_top_3D = repeat(reshape(depth_top, (1,1,ndepth)), outer=(nlat,nlon,1))
A_2D = permutedims(vars["area2x2"], (2,1))m^2
nboxes = Int(vars["nOCEAN"])
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


println("  Circulation")
iwet_orig = findall(vec(wet3D_orig))
permuted_iwet = permutedims(LinearIndices(size(wet3D)), (2,1,3))[iwet_orig]
p = sortperm(permuted_iwet)
#test that sort `p` is correct
n = prod(collect(size(wet3D)))
x3D_orig = wet3D_orig .* reshape(1:n, size(wet3D_orig))
x_orig = x3D_orig[iwet_orig]
# permute x
x = x_orig[p]
# now reconstruct x3d_orig from x
x3D = Int.(wet3D)
x3D[iwet] .= x
println(permutedims(x3D, (2,1,3)) == x3D_orig)
v = ustrip.(volume_3D[iwet])
T = (-1/s) * sparse(Diagonal(inv.(v))) * (vars["A_adv"] + vars["A_dif"])[p,p] # `p` swaps lat and lon? TODO check

name = "OCCA"
data_path = "/Users/benoitpasquier/Data"
bson_dir = joinpath(data_path, "OceanGrids")
println("Saving as BSON file in $bson_dir")
bson_file = joinpath(bson_dir, "$name.bson")
isdir(bson_dir) || mkdir(bson_dir)
isfile(bson_file) && rm(bson_file)

BSON.@save bson_file grid T


