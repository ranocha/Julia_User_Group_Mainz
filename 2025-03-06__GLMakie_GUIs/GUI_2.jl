#GUI_2
using GLMakie

a = Observable(1.0)
b = Observable(1.0)

f1(x) = a[].*x.^3 + b[].*x.^2
f2(x) = a[].*x.^2 + b[].*x

fig = Figure()

menu = Menu(fig, options = ["red", "blue", "green","brown"], default = "red")

funcs = [f1, f2]

menu2 = Menu(fig,
    options = zip(["ax³ + bx²","ax² + bx"], funcs),
    default = "ax³ + bx²")

text_a = Textbox(fig, placeholder = "1.0", validator = Float64)
text_b = Textbox(fig, placeholder = "1.0", validator = Float64)

fig[1, 1] = vgrid!(
    Label(fig, "Colormap", width = nothing),
    menu,
    Label(fig, "Function", width = nothing),
    menu2,
    Label(fig, "a:", width = 200),
    text_a,
    Label(fig, "b:", width = 200),
    text_b;
    tellheight = false, width = 200)

ax = Axis(fig[1, 2])

func = Observable{Any}(funcs[1])

xs = -2:0.01:2
ys = lift(func) do f
    f.(xs)
end
line = lines!(ax, xs, ys, color = :red)

on(menu.selection) do s
    line.color = s
end
notify(menu.selection)


on(menu2.selection) do s
    func[] = s
    autolimits!(ax)
end
notify(menu2.selection)

on(text_a.stored_string) do s
    a[] = parse(Float64, s)
    notify(a)
    notify(menu2.selection)
end

on(text_b.stored_string) do s
    b[] = parse(Float64, s)
    notify(b)
    notify(menu2.selection)
end

fig