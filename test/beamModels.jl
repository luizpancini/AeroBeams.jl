using AeroBeams, LinearAlgebra

# simpleBeamModel
beam1 = Beam(name="beam1",length=2.0,nElements=5,normalizedNodalPositions=[0.0,0.1,0.25,0.45,0.7,1.0],f_A_of_x1t=(x1,t)->[1.0*x1*t,1,2.0*x1],C=[Matrix(1.0*LinearAlgebra.I,6,6)])

pointMass = PointInertia(elementLocalID=3,mass=1.0,Ixx=1.0)

add_point_inertias_to_beam!(beam1,inputPointInertias=[pointMass])

# update_beam!(beam=beam1)

clamp1 = create_BC(name="clamp1",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

forces1 = create_BC(name="tipForces",beam=beam1,node=6,types=["F1A","F3A"],values=[10.0,(t)->1.0-t])

simpleBeamModel = Model(name="simpleBeam",beams=[beam1],BCs=[clamp1])

# splitBeamModel
beam2 = Beam(name="beam2",connectedBeams=[beam1],connectedNodesThis=[1],connectedNodesOther=[5],length=0.5,nElements=4,C=[Matrix(1.0*LinearAlgebra.I,6,6)],rotationParametrization="E321",p0=[0.0;-π/2;0.0],u0_of_x1=x1->[0.0*x1/0.5;0.0;0.0])

clamp2 = create_BC(name="clamp2",beam=beam2,node=5,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

add_initial_displacements_and_velocities_to_beam!(beam2,conditionTypes=["p0_of_x1"],conditionFuns=[x1->[0.01*x1/0.5;0.0;0.0]])

add_loads_to_beam!(beam2,loadTypes=["m_A_of_x1t"],loadFuns=[(x1,t)->[1.0*x1*t,1,2.0*x1]])

splitBeamModel = Model(name="splitBeam")
splitBeamModel.beams = [beam1,beam2]
splitBeamModel.BCs = [clamp1,clamp2]
update_model!(splitBeamModel)

# splitBeamModel 2
splitBeamModel2 = Model(name="splitBeam2",beams=[beam1,beam2],BCs=[clamp1,clamp2])

plot_undeformed_assembly(splitBeamModel2)

# Create and solve problem
# problem = SteadyProblem(model=splitBeamModel)
# solve!(problem)

println("Finished")
