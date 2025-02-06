using GLMakie
using HOHQMesh

sample = newProject("sample","out")

addBackgroundGrid!(sample, [2.0, 2.0, 0.0])
setMeshFileFormat!(sample, "ABAQUS")
setPolynomialOrder!(sample ,3)
radius = 5.0

L = 20.0

inflow = newCircularArcCurve("inflow",
			     [0.0, 0.0, 0.0], # center
			     radius,		      # radius
			     90.0,	      # start angle	
			     180.0,	      # end angle
			     "degrees")

cylinder = newCircularArcCurve("cylinder",
			       [0.5, 0.0, 0.0],
			       0.5,
			       180.0,
			       0.0,
			       "degrees")

outflow = newEndPointsLineCurve("outflow", [L, 0.0, 0.0], [L, radius, 0.0])

top = newEndPointsLineCurve("top", [L, radius, 0.0], [0.0, radius, 0.0])

symmetry_line1 = newEndPointsLineCurve(":symmetry", [-radius, 0.0, 0.0], [0.0, 0.0, 0.0])

symmetry_line2 = newEndPointsLineCurve(":symmetry", [1.0, 0.0, 0.0], [L, 0.0, 0.0])

addCurveToOuterBoundary!(sample, inflow)
addCurveToOuterBoundary!(sample, symmetry_line1)
addCurveToOuterBoundary!(sample, cylinder)
addCurveToOuterBoundary!(sample, symmetry_line2)
addCurveToOuterBoundary!(sample, outflow)
addCurveToOuterBoundary!(sample, top)

plotProject!(sample, MODEL+GRID)

wake_region = newRefinementLine("region","smooth",[1.0, 0.0, 0.0] , [L, 0.0, 0.0], 0.2, 1.0)

addRefinementRegion!(sample, wake_region)

cylinder_ref = newRefinementCenter("region", "smooth", [0.5, 0.0, 0.0], 0.2, 1.0)

addRefinementRegion!(sample, cylinder_ref)

generate_mesh(sample)




