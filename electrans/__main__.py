import argparse
import os
from . import transition_analysis_utils as tau
from . import vtk_output_utils as vtk_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir",
                        type=str,
                        help="The directory containing the input cube files for transitions. "
                        "The directory should a 'metadata.csv' file which contains 3 columns. "
                        "The first specifies the name of the transition. The second and third "
                        "columns specify hole and particle cube files respectively.")
    parser.add_argument("-p", "--different_atom_pos",
                        default=False,
                        action='store_true',
                        help="If specified, Voronoi diagram will be computed for each cube "
                        "file separately. By default it is assumed, the atomic positions remain "
                        "the same for all the cube files.")
    parser.add_argument("-q", "--use_quadratic_optimization",
                        default=False,
                        action='store_true',
                        help="Use quadratic programming based optimization approach to compute "
                        "the charge trasfer between subgroups. By default a hueristic approach "
                        "is used for this purpose.")
    parser.add_argument("-d", "--output_transition_diagram",
                        default=False,
                        action='store_true',
                        help="If specified, transition diagram will be generated for all the "
                        "transitions in '<input_dir>/results/transition_diagrams/' directory.")
    parser.add_argument("-a", "--output_atomic_charges",
                        default=False,
                        action='store_true',
                        help="If specified, atomic charges for each transition will be output"
                        " as a csv file in the directory '<input_dir>/results/atomic_charges/'.")
    parser.add_argument("-s", "--output_subgroup_charges",
                        default=False,
                        action='store_true',
                        help="If this option is set, subgroup charges and amount of charge "
                        "transfer for each transition would be output in a directory named "
                        "'<input_dir>/results/subgroup_charges/'.")
    parser.add_argument("-av", "--output_atoms_vtk",
                        default=False,
                        action='store_true',
                        help="If this option is set, the atoms are saved as 3D model in VTK format. "
                        "A file is generated for each transition which contains the data about "
                        "the hole and particle charge, charge difference, subgroup, etc. "
                        "These VTK files are saved in '<input_dir>/results/vtk/atoms/' and "
                        "can be loaded in VTK compatible software like Paraview.")
    parser.add_argument("-sv", "--output_segmentation_vtk",
                        default=False,
                        action='store_true',
                        help="If specified, the computed Voronoi segmentation at atomic and "
                        "subgroup scales are saved as VTK compatible files in the directory "
                        "'<input_dir>/results/vtk/segmentation/'. These files can be loaded in "
                        "VTK compatible software like Paraview.")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=4,
                        help="Specify the number of parallel threads to be used in "
                        "computing the Voronoi diagram. By default four threads are used.")

    args = parser.parse_args()
    input_dir = args.input_dir
    different_atom_pos = args.different_atom_pos
    use_quadratic_optimization = args.use_quadratic_optimization
    threads = args.threads
    output_transition_diagram = args.output_transition_diagram
    output_atomic_charges = args.output_atomic_charges
    output_subgroup_charges = args.output_subgroup_charges
    output_atoms_vtk = args.output_atoms_vtk
    output_segmentation_vtk = args.output_segmentation_vtk

    if not os.path.isdir(input_dir):
        print("Can't find the specified input directory. Exiting ...")
        exit(0)

    metadata_file = input_dir + "metadata.csv"
    if not os.path.isfile(metadata_file):
        print("Can't find metadata.csv file in the directory. Exiting ...")
        exit(0)

    state_files = tau.read_metadata_file(metadata_file)
    if not state_files:
        print("No cube files specified in the metadata.csv file. Exiting ...")
        exit(0)

    output_dir = input_dir + "results/"
    if not os.path.isdir(output_dir):
        print("Creating the output directory: %s" % output_dir)
        os.makedirs(os.path.dirname(output_dir))

    if output_atomic_charges:
        output_dir = input_dir + "results/atomic_charges/"
        if not os.path.isdir(output_dir):
            print("Creating the output directory for atomic charges: %s" %
                  output_dir)
            os.makedirs(os.path.dirname(output_dir))

    if output_subgroup_charges:
        output_dir = input_dir + "results/subgroup_charges/"
        if not os.path.isdir(output_dir):
            print("Creating the output directory for subgroup charges: %s" %
                  output_dir)
            os.makedirs(os.path.dirname(output_dir))

    if output_transition_diagram:
        output_dir = input_dir + "results/transition_diagrams/"
        if not os.path.isdir(output_dir):
            print("Creating the output directory for transition diagrams: %s" %
                  output_dir)
            os.makedirs(os.path.dirname(output_dir))

    if output_atoms_vtk:
        output_dir = input_dir + "results/vtk/atoms/"
        if not os.path.isdir(output_dir):
            print("Creating the output directory for atomic charges: %s" %
                  output_dir)
            os.makedirs(os.path.dirname(output_dir))

    transitions = []
    transition_names = []
    for state_file in state_files:
        hole_cubeFile = input_dir + state_file[1]
        particle_cubeFile = input_dir + state_file[2]
        transition = tau.load_transition(hole_cubeFile, particle_cubeFile)
        transitions.append(transition)
        transition_names.append(state_file[0])

    print("Loaded the hole and particle cube files ...")

    segment_arrays = tau.compute_atomic_charges(
        transitions, num_threads=threads, same_atomic_positions=not different_atom_pos,
        save_segmention=output_segmentation_vtk)
    print("Compted the Voronoi diagram based segmentation and atomic charges ...")

    if output_subgroup_charges or output_transition_diagram:
        subgroup_file = input_dir + "subgroups.txt"
        if os.path.isfile(metadata_file):
            subgroup_names, atom_subgroup_map = tau.load_subgroups(
                subgroup_file)
        else:
            print("Can't find subgroups.txt file in the input directory. "
                  "The subgroup charges and transition diagrams will not be computed.")
            output_subgroup_charges = output_transition_diagram = False
            subgroup_names = []
            atom_subgroup_map = [0] * transitions[0].num_atoms()

    for i in range(len(transitions)):
        transition = transitions[i]
        if output_atomic_charges:
            output_file = input_dir + "results/atomic_charges/" + \
                "%s.csv" % transition_names[i]
            tau.save_atomic_charges(
                output_file, transition.hole_data.atoms, transition.hole_charges,
                transition.particle_charges, subgroup_names, atom_subgroup_map)

        if output_atoms_vtk:
            output_file = input_dir + "results/vtk/atoms/" + \
                "%s.vtp" % transition_names[i]
            vtk_out.write_atoms(output_file, transition.hole_data.atoms,
                                transition.hole_charges, transition.particle_charges, atom_subgroup_map)

        if output_subgroup_charges or output_transition_diagram:
            subgroup_info = tau.SubgroupInfo()
            subgroup_info.set_subgroups(subgroup_names, atom_subgroup_map)
            tau.compute_subgroup_charges(transition, subgroup_info, use_hueristic=not use_quadratic_optimization)
            if output_subgroup_charges:
                output_file = input_dir + "results/subgroup_charges/" + \
                    "%s.txt" % transition_names[i]
                subgroup_info.save_to_file(output_file)
            if output_transition_diagram:
                output_file = input_dir + "results/transition_diagrams/" + \
                    "%s.pdf" % transition_names[i]
                tau.create_diagram(
                    subgroup_info, title=transition_names[i], show_plot=False, save_plot=True, file_name=output_file)

    if output_segmentation_vtk:
        output_dir = input_dir + "results/vtk/segmentation/"
        if not os.path.isdir(output_dir):
            os.makedirs(os.path.dirname(output_dir))
        if different_atom_pos:
            for i in range(len(transitions)):
                data = transitions[i].hole_data
                vtk_out.write_segments(output_dir + "%s_hole_seg_atoms.vtp" %
                                       transition_names[i], segment_arrays[2 * i], data.basis, data.atoms)
                if output_subgroup_charges or output_transition_diagram:
                    vtk_out.write_subgroup_segments(
                        output_dir +
                        "%s_hole_seg_subgroups.vtp" % transition_names[i],
                        segment_arrays[2 * i], data.basis, data.atoms, len(subgroup_names), atom_subgroup_map)

                data = transitions[i].particle_data
                vtk_out.write_segments(output_dir + "%s_particle_seg_atoms.vtp" %
                                       transition_names[i], segment_arrays[2 * i + 1], data.basis, data.atoms)
                if output_subgroup_charges or output_transition_diagram:
                    vtk_out.write_subgroup_segments(
                        output_dir +
                        "%s_particle_seg_subgroups.vtp" % transition_names[i],
                        segment_arrays[2 * i + 1], data.basis, data.atoms, len(subgroup_names), atom_subgroup_map)
        else:
            data = transitions[0].hole_data
            vtk_out.write_segments(
                output_dir + "seg_atoms.vtp", segment_arrays[0], data.basis, data.atoms)
            if output_subgroup_charges or output_transition_diagram:
                vtk_out.write_subgroup_segments(
                    output_dir + "seg_subgroups.vtp", segment_arrays[0], data.basis, data.atoms, len(subgroup_names), atom_subgroup_map)


if __name__ == "__main__":
    main()
