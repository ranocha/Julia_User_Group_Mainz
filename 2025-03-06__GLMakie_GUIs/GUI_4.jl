#GUI_4 - draw a curve by clicking repeatedly on the axis. Press enter to display a list of points
using GLMakie

fig = Figure()
ax = Axis(fig[1,1])

for interact in keys(ax.interactions)
    deactivate_interaction!(ax, interact)
end
points = Observable(Point2f[])
lines!(ax, points, color=:black)
scatter!(ax, points, color=:gray, markersize=20) 

on(events(ax.scene).mousebutton) do event
    if event.button == Mouse.left
      if event.action == Mouse.press 
        mp = mouseposition(ax.scene)
        push!(points[], mp)
        notify(points)
      end
    end
end

on(events(ax.scene).keyboardbutton) do event
  if event.key == Keyboard.enter
    println("pressed enter")
    @show points[]
  end
end

display(fig)