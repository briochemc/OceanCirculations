#=
This is to test the ideal age computation
=#

using SparseArrays          # For sparse matrix
using MAT                   # For loading MAT format
using LinearAlgebra

# No download, this is a local test file

println("Reading grid file (grid data)")
data_path = ""
mat_boxes_file = joinpath(data_path, "MITgcm_2.8deg/Matrix5/Data/boxes.mat")
boxes = matread(mat_boxes_file)
mat_grd_file = joinpath(data_path, "MITgcm_2.8deg/grid.mat")
grd = matread(mat_grd_file)

println("Reading Matrix file (circulation matrix)")
mat_Matrix_file = joinpath(data_path, "MITgcm_2.8deg/Matrix5/TMs/matrix_nocorrection_annualmean.mat")
vars = matread(mat_Matrix_file)

println("  Wet boxes")
wet3D = convert(BitArray{3}, grd["bathy"])
iwet = findall(vec(wet3D))

println("  Circulation (T)")
Te = vars["Aexpms"]
Ti = vars["Aimpms"]

#=
The equation for the ideal mean age is

∂a/∂t + T a = 1 - Λ a

where Λ = Diagonal(z == z_surface) / τ and τ = 1s

But we don't have T here, using the TMM.
With the TMM, we have Ti and Te, and the time stepping is done via

a(i+1) = Ti (Te a(i) + τ_T sms)

So the steady-state solution is given by

(I - Ti Te) a = τ_T Ti sms

Here sms = 1 - Λ a, so we solve

(I - Ti Te + τ_t Ti Λ) a = τ_T Ti 1

or 

a = M \ rhs
=#
τ_T = grd["deltaT"]
z = vec(boxes["Zbox"])
Λ = sparse(Diagonal(z .< 30)) / 1

M = I - Ti * Te + τ_T * Ti * Λ

nb = length(z)
rhs = τ_T * Ti * ones(nb)

a = M \ rhs

#=
With operator splitting
=#

i = findall(z .> 30) # indices of interior boxes

a = (Ti[i,i] * Te[i,i] - I) \ (-τ_T * Ti[i,i] * ones(length(i)))



println("  Grid")
lat = vec(grd["y"]) * u"°"
lon = vec(grd["x"]) * u"°"
depth = vec(grd["z"]) * u"m"
nlat, nlon, ndepth = Int(grd["ny"]), Int(grd["nx"]), Int(grd["nz"])
δx_3D = repeat(reshape(permutedims(grd["dx"], [2,1]), (nlat,nlon,1)), outer=(1,1,ndepth)) * u"m"
δy_3D = repeat(reshape(permutedims(grd["dy"], [2,1]), (nlat,nlon,1)), outer=(1,1,ndepth)) * u"m"
δz_3D = repeat(reshape(grd["dznom"], (1,1,ndepth)), outer=(nlat,nlon,1)) * u"m"
A_2D = permutedims(grd["da"], [2,1,3])[:,:,1] * u"m^2"   # only on iwet!
# TODO Figure out how the volume is computed!
volume_3D = permutedims(grd["da"] .* grd["dz"], [2,1,3]) # only on iwet!
volume_3D = δx_3D .* δy_3D .* δz_3D

# Also the direction? (Divergence or convergence)?

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



bson_dir = joinpath(data_path, "OceanGrids")
println("Saving as BSON file in $bson_dir")
bson_file = joinpath(bson_dir, "POPX3.bson")
isdir(bson_dir) || mkdir(bson_dir)
isfile(bson_file) && rm(bson_file)
BSON.@save bson_file grid wet3D T
