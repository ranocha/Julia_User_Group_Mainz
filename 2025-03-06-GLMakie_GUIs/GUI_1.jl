# First GUI


# 1) Simple bar plot
x = [1,2,3]
y = [2,1,3]

using GLMakie
fig = Figure()
ax = Axis(fig[1,1])
barplot!(ax, x, y)
display(fig)


# 2) Same plot with y as Observable:
x = [1,2,3]
y = Observable([2,1,3])

using GLMakie
fig = Figure()
ax = Axis(fig[1,1])
barplot!(ax, x, y)
display(fig)

y[] = [1,5,4] # update y
#display(fig)  # show figure again  

# 3) add a button
x = [1,2,3]
y = Observable([2,1,3])

using GLMakie
fig = Figure()
ax  = Axis(fig[1,1:3])
but = Button(fig[2,2], label="y[2]") 
barplot!(ax, x, y)

display(fig)

# 4) Give the button some action
x = [1,2,3]
y = Observable([2,1,3])

using GLMakie
fig = Figure()
ax  = Axis(fig[1,1:3])
but = Button(fig[2,2], label="y[2]") 
barplot!(ax, x, y)

on(but.clicks) do n
    @show n
    y[][2] +=1    # update y[2]
    notify(y)     # update plot
end

on(y) do yval
    @show yval
    println("y changed")
end
display(fig)
