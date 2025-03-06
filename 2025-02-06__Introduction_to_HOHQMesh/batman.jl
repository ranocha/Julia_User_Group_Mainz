using GLMakie, HOHQMesh

batman_project = newProject("BatmanLogo", "out")
addBackgroundGrid!(batman_project, [0.6, 0.6, 0.0])
setPolynomialOrder!(batman_project, 3)

function spline_batman(t)
    x =
        -0.5 * t * 21 + 0.25 * abs(t * 21 - 1) - 0.5 * abs(t * 21 - 3) + 0.75 * abs(t * 21 - 5) + 1.5 * abs(t * 21 - 13) - 0.25 * abs(t * 21 - 17) - 1.25 * abs(t * 21 - 21) - 7 * sin((pi / 12) * (abs(t * 21 - 5) - abs(t * 21 - 8) + 3)) -
        (1 / 100) * (abs(t * 21 - 8) - abs(t * 21 - 13) - 5)^3 - 1.5
    y =
        3 / 4 * abs(t * 21 - 1) - (3 / 4) * abs(t * 21 - 3) - (7 / 5) * abs(t * 21 - 8) + (7 / 5) * abs(t * 21 - 13) + (7 / 16) * (abs(t * 21 - 3) - abs(t * 21 - 5) - 2)^2 + 4 * sin((pi / 12) * (abs(t * 21 - 5) - abs(t * 21 - 8) - 3)) -
        (5 / 16) * (abs(t * 21 - 13) - abs(t * 21 - 17))^2 - (1 / 4) * (abs(t * 21 - 17) - abs(t * 21 - 21) + 2)^2 + 11.5
    return [t x y 0]
end

t = range(0,1,1001)
spline_batman_data = zeros(round(length(t)), 4)
for i in eachindex(t)
    spline_batman_data[i, :] = spline_batman(t[i])
end

batman_spline = newSplineCurve("Batman", round(length(t)), spline_batman_data)
symmetry_line = newEndPointsLineCurve(":symmetry", [0.0, -8.0, 0.0], [0.0, 14.0, 0.0])

addCurveToOuterBoundary!(batman_project, symmetry_line)
addCurveToOuterBoundary!(batman_project, batman_spline)

plotProject!(batman_project, MODEL + GRID)

generate_mesh(batman_project)