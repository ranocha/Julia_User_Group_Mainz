using HOHQMesh
using GLMakie
using HTTP

function load_airfoil(airfoil_name::String)
    url = "http://airfoiltools.com/airfoil/seligdatfile?airfoil=$(airfoil_name)"
    
    response = HTTP.get(url)

    righe = split(String(response.body), '\n')

    airfoil = []

    for riga in righe
        colonne = split(riga)
        if length(colonne) == 2
            push!(airfoil, (parse(Float64, colonne[1]), parse(Float64, colonne[2])))
        end
    end

    airfoil_coordinates = hcat([x[1] for x in airfoil], [y[2] for y in airfoil])

    x = airfoil_coordinates[:, 1]
    y = airfoil_coordinates[:, 2]

    if x[1] != x[end] || y[1] != y[end]
        x = vcat(x, x[1])
        y = vcat(y, y[1])
    end

    t = [0.0]
    for i in 2:length(x)
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        distance = sqrt(dx^2 + dy^2)
        push!(t, t[end] + distance)
    end
    
    t .= t ./ t[end]
    
    result = hcat(t, x, y, zeros(length(x)))
    return result

end

airfoil_project = newProject("airfoil", "out")
addBackgroundGrid!(airfoil_project, [1.0, 1.0, 0.0])
setPolynomialOrder!(airfoil_project, 4)
L = 2.0
inflow = newCircularArcCurve("inflow",       # curve name
	[0.0, 0.0, 0.0], # circle center
	L,             # circle radius
	90.0,             # start angle
	270.0,           # end angle
	"degrees")

bottom = newEndPointsLineCurve("Bottom",          # curve name
	[0.0, -L, 0.0], # start point   
	[L, -L, 0.0]) # end point                              

top = newEndPointsLineCurve("Top",          # curve name
	[L, L, 0.0], # start point
	[0.0,  L, 0.0]) # end point 

outflow = newEndPointsLineCurve("Outflow",          # curve name
	[L, -L, 0.0], # start point
	[L, L, 0.0]) # end point 

addCurveToOuterBoundary!(airfoil_project, inflow)
addCurveToOuterBoundary!(airfoil_project, bottom)
addCurveToOuterBoundary!(airfoil_project, outflow)
addCurveToOuterBoundary!(airfoil_project, top)

airfoil_name = "b737c-il"
airfoil_spline = newSplineCurve("Airfoil", size(load_airfoil(airfoil_name))[1], load_airfoil(airfoil_name))

addCurveToInnerBoundary!(airfoil_project, airfoil_spline, "Inner")

plotProject!(airfoil_project, MODEL + GRID)

# Refinements
leading_edge = newRefinementCenter("region", "smooth", [0.2, 0.0, 0.0], 0.01, 0.3)
trailing_edge = newRefinementCenter("region", "smooth", [0.8, 0.0, 0.0], 0.01, 0.3)
center = newRefinementCenter("region", "smooth", [0.5, 0.0, 0.0], 0.01, 0.3)

addRefinementRegion!(airfoil_project, leading_edge)
addRefinementRegion!(airfoil_project, trailing_edge)
addRefinementRegion!(airfoil_project, center)

generate_mesh(airfoil_project)