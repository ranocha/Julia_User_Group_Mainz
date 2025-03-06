using HOHQMesh
using GLMakie

tutorial = newProject("box_two_circles", "out")

setPolynomialOrder!(tutorial, 4)
setMeshFileFormat!(tutorial, "ABAQUS")

bounds = [15.0, 0.0, 0.0, 30.0]
N = [30, 15, 0]
addBackgroundGrid!(tutorial, bounds, N)

circle1 = newCircularArcCurve("circle1",       # curve name
                              [4.0, 4.0, 0.0], # circle center
                              2.0,             # circle radius
                              0.0,             # start angle
                              360.0,           # end angle
                              "degrees") 

circle2 = newCircularArcCurve("circle2",        # curve name
                              [20.0, 9.0, 0.0], # circle center
                              4.0,              # circle radius
                              0.0,              # start angle
                              2.0 * pi,         # end angle
                              "radians")        # angle units

addCurveToInnerBoundary!(tutorial, circle1, "inner1")                              
addCurveToInnerBoundary!(tutorial, circle2, "inner2")

plotProject!(tutorial, MODEL + GRID)

generate_mesh(tutorial)