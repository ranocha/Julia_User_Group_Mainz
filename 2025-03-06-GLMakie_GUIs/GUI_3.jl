#GUI_3 - feedback after clicking on the axis
using GLMakie

fig = Figure()
ax = Axis(fig[1,1])

deactivate_interaction!(ax, :scrollzoom)
deactivate_interaction!(ax, :limitreset)
deactivate_interaction!(ax, :dragpan)
deactivate_interaction!(ax, :rectanglezoom)

vlines!(ax, 0.5, linewidth=10)
hlines!(ax, 0.5, linewidth=10)

point = Observable(Point2(0.5,0.5))
scatter!(ax, point, color=:red, markersize=10) 

on(events(ax.scene).mousebutton) do e
    position = mouseposition(ax.scene)
    point[] = position
    @show position
end


display(fig)