# GUI to plot Topography
using GeophysicalModelGenerator, GMT, GLMakie

# Download 

# Example of loading Topo:
Topo = import_topo([4,20,37,49]);


# Plot Topo with heatmap(Topo) or with:
heatmap(Topo.lon.val[:,1,1],Topo.lat.val[1,:,1], ustrip.(Topo.fields.Topography[:,:,]))


# Exercise: create a GUI around this...
