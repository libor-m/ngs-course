# the idea is to present the basic idea of pca with a nice example
# taking a 3d banana, rotaing it so it's indistinguishable from 
# the main axes, and then letting pca to find the best rotation

# get a banana at 
# http://www.turbosquid.com/FullPreview/Index.cfm/ID/458771

# .obj maybe can be converted to a point cloud pretty easily
<Banana.obj grep '^v ' | cut -d' ' -f2- > banana.tsv
# http://www.mlahanas.de/CompGeom/Point_in_Poly.htm

# TODO: make a shiny app?
# upload to https://www.shinyapps.io/?