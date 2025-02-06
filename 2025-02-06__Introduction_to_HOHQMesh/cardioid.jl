using GLMakie, HOHQMesh

cardioid_project = newProject("TheBlob", "out")
a = 8
xEqn = "x(t) = 2 * $a * ( 1 - cos(2 * pi * t))* cos(2 * pi * t)"
yEqn = "y(t) = 2 * $a * ( 1 - cos(2 * pi * t))* sin(2 * pi * t)"
zEqn = "z(t) = 0.0"
cardioid = newParametricEquationCurve("Blob", xEqn, yEqn, zEqn)

addCurveToOuterBoundary!(cardioid_project, cardioid)

addBackgroundGrid!(cardioid_project, [2.0, 2.0, 0.0])

plotProject!(cardioid_project, MODEL+GRID)

setBackgroundGridSize!(cardioid_project, 0.5, 0.5)

generate_mesh(cardioid_project)