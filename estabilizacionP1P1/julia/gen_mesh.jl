using Plots
using Gmsh

gmsh.initialize()

#…define the model

model = gmsh.model

model.add("Square")

#…define mesh density near a point

cl = 0.1;

#…define four points in the geometry

gmsh.model.geo.addPoint(0.5, 0.5, 0., cl, 1)

gmsh.model.geo.addPoint(-0.5, 0.5, 0., cl, 2)

gmsh.model.geo.addPoint(-0.5, -0.5,0., cl, 3)

gmsh.model.geo.addPoint(0.5, -0.5, 0., cl, 4)

#…define four edges in the geometry

gmsh.model.geo.addLine(1, 2, 1)

gmsh.model.geo.addLine(2, 3, 2)

gmsh.model.geo.addLine(3, 4, 3)

gmsh.model.geo.addLine(4, 1, 4)

#…define outer boundary

gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

#…define planar surface

gmsh.model.geo.addPlaneSurface([1], 1)

#…define physics

gmsh.model.addPhysicalGroup(1, [1], 1)

#…synchronize the model

gmsh.model.geo.synchronize()

#…generate the mesh in 2D

model.mesh.generate(2)

#…save the mesh to file for future reference

gmsh.write("square.msh")

#######

#gmsh.finalize()

node_ids, node_coord, _ = gmsh.model.mesh.getNodes()

nnodes = length(node_ids)

xnode = node_coord[1:3:end]

ynode = node_coord[2:3:end]

#…Plot the mesh nodes

#…plots nodes only

#scatter(xnode, ynode)


#…or alternatively

using GR

z = ones(length(xnode))

trisurf(xnode,ynode,z)


#…Retrieve the set of mesh nodes and the number of nodes (nnodes)

#…sort the node coordinates by ID, such that Node one sits at row 1

tosort = [node_ids node_coord[1:3:end] node_coord[2:3:end]];

sorted = sortslices(tosort , dims = 1);

node_ids = sorted[:,1]

xnode = sorted[:,2]

ynode = sorted[:,3]