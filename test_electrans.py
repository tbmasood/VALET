from electrans import transition_analysis_utils as tau

def test_example():
    hole_cubeFile = "./data/tq-td_State1-GS.cube"
    particle_cubeFile = "./data/tq-td_State1-ES.cube"

    # Load the the two cube files corresponding to hole (ground state) and particle (excited state)
    transition = tau.load_transition(hole_cubeFile, particle_cubeFile)
    print("Loaded the hole and particle cube files ...")

    # Compute the weighted Voronoi diagram and accumulate charges for each atom
    segment_array = tau.compute_atomic_charges([transition])
    print("Compted the Voronoi diagram based segmentation and atomic charges ...")
    
    # Provide the grouping of atoms in meaningful subgroups
    subgroup_names = ["THIO", "QUIN"]
    atom_subgroup_map = [-1] * transition.num_atoms()
    # In this case the first 8 atoms belong to Thiophene while the rest of the atoms belong to Quinoxaline
    atom_subgroup_map[:8] = [0] * 8
    atom_subgroup_map[8:] = [1] * (transition.num_atoms() - 8)
    subgroup_info = tau.SubgroupInfo()
    subgroup_info.set_subgroups(subgroup_names, atom_subgroup_map)
    
    # Compute the charges and the amount of charge transfer between subgroups
    tau.compute_subgroup_charges(transition, subgroup_info)
    print("Compted subgroup charges and charge trasfer ...")
    
    # Generate and plot the Transition Diagram
    tau.create_diagram(subgroup_info, title="TQ (State 1)")

    # Output VTK files
    from electrans import vtk_output_utils as vtk_out

    data = transition.hole_data
    vtk_out.write_segmentation("dens_0_GS.vti", data.scalars, data.basis)
    vtk_out.write_segmentation("dens_0_ES.vti", transition.particle_data.scalars, data.basis)
    vtk_out.write_atoms("atoms_0_GS.vtp", data.atoms, transition.hole_charges)
    vtk_out.write_atoms("atoms_0_ES.vtp", data.atoms, transition.particle_charges)
    chargeDiff = [0] * data.num_atoms()
    for atomID in range(data.num_atoms()):
        chargeDiff[atomID] = transition.particle_charges[atomID] - transition.hole_charges[atomID]
    vtk_out.write_atoms("atoms_0_diff.vtp", data.atoms, chargeDiff)
    vtk_out.write_segmentation("seg.vti", segment_array, data.basis, array_name="segment")
    vtk_out.write_segments("seg.vtp", segment_array, data.basis, data.atoms)
    vtk_out.write_subgroup_segments("seg_ligands.vtp", segment_array, data.basis, data.atoms, [range(8), range(8, data.num_atoms())])

    print("All output VTK files written to disk ...")

if __name__ == "__main__":
    test_example()