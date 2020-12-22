from electrans import transition_analysis_utils as tau


def test_example():
    hole_cubeFile = "./data/tq-td_State1-GS.cube"
    particle_cubeFile = "./data/tq-td_State1-ES.cube"

    # Load the the two cube files corresponding to hole (ground state) and particle (excited state)
    transition = tau.load_transition(hole_cubeFile, particle_cubeFile)
    print("Loaded the hole and particle cube files ...")

    # Compute the weighted Voronoi diagram and accumulate charges for each atom
    segment_arrays = tau.compute_atomic_charges([transition], save_segmention=True)
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
    vtk_out.write_segmentation(
        "dens_0_ES.vti", transition.particle_data.scalars, data.basis)
    vtk_out.write_atoms("atoms.vtp", data.atoms, transition.hole_charges,
                        transition.particle_charges, atom_subgroup_map)
    vtk_out.write_segmentation(
        "seg.vti", segment_arrays[0], data.basis, array_name="segment")
    vtk_out.write_segments(
        "seg.vtp", segment_arrays[0], data.basis, data.atoms)
    vtk_out.write_subgroup_segments("seg_ligands.vtp", segment_arrays[0], data.basis, data.atoms, len(
        subgroup_names), atom_subgroup_map)

    print("All output VTK files written to disk ...")


if __name__ == "__main__":
    test_example()
