class ElectronicTransition:
    def __init__(self, hole_data, particle_data):
        self.hole_data = hole_data
        self.particle_data = particle_data
        self.hole_charges = [0] * hole_data.num_atoms()
        self.particle_charges = [0] * particle_data.num_atoms()
        if hole_data.num_atoms() != particle_data.num_atoms():
            raise Exception(
                "The number of atoms in Hole and particle cube files are different")

    def num_atoms(self):
        return self.hole_data.num_atoms()

    def normalize(self):
        sumHoleCharges = sum(self.hole_charges)
        sumParticleCharges = sum(self.particle_charges)
        for i in range(self.num_atoms()):
            self.hole_charges[i] = 2.0 * self.hole_charges[i] / sumHoleCharges
            self.particle_charges[i] = 2.0 * \
                self.particle_charges[i] / sumParticleCharges


def load_transition(hole_cubefile: str, particle_cubefile: str):
    from . import cubefile_utils as cubeFileUtils

    hole_data = cubeFileUtils.load_cubefile(hole_cubefile)
    particle_data = cubeFileUtils.load_cubefile(particle_cubefile)

    return ElectronicTransition(hole_data, particle_data)


def _voronoi_subprocess(x, basis, atoms, size):
    import numpy

    segment_array = numpy.empty((size[1], size[2]), dtype=numpy.int32)
    for y in range(size[1]):
        for z in range(size[2]):
            pt = (x * basis[0][0], y * basis[1][1], z * basis[2][2])
            closest = 0
            minDist = float('inf')
            for i in range(len(atoms)):
                xDiff = atoms[i].coordinates[0] - pt[0]
                yDiff = atoms[i].coordinates[1] - pt[1]
                zDiff = atoms[i].coordinates[2] - pt[2]
                dist = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff \
                    - (atoms[i].radius * atoms[i].radius)
                if dist < minDist:
                    closest = i
                    minDist = dist
            segment_array[y][z] = closest
    return segment_array


def _voronoi_parallel(basis, atoms, size, num_threads):
    import numpy
    from multiprocessing import Pool
    from functools import partial

    segment_array = numpy.empty(size, dtype=numpy.int32)
    with Pool(num_threads) as pool:
        ret = pool.map(partial(_voronoi_subprocess,
                               basis=basis, atoms=atoms, size=size), range(size[0]))
        for x in range(size[0]):
            segment_array[x] = ret[x]
    return segment_array


def _voronoi_sequential(basis, atoms, size):
    import numpy

    segment_array = numpy.empty(size, dtype=numpy.int32)
    for x in range(size[0]):
        for y in range(size[1]):
            for z in range(size[2]):
                pt = (x * basis[0][0], y * basis[1][1], z * basis[2][2])
                closest = 0
                minDist = float('inf')
                for i in range(len(atoms)):
                    xDiff = atoms[i].coordinates[0] - pt[0]
                    yDiff = atoms[i].coordinates[1] - pt[1]
                    zDiff = atoms[i].coordinates[2] - pt[2]
                    dist = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff \
                        - (atoms[i].radius * atoms[i].radius)
                    if dist < minDist:
                        closest = i
                        minDist = dist
                segment_array[x][y][z] = closest
    return segment_array


def voronoi_segmentation(basis, atoms, size, num_threads):
    if num_threads < 2:
        return _voronoi_sequential(basis, atoms, size)
    else:
        return _voronoi_parallel(basis, atoms, size, num_threads)


def accumulate_atomic_charges(segment_array, scalars, atoms, size):
    chargePerAtom = [0] * len(atoms)
    atomVolume = [0] * len(atoms)
    totalSum = 0
    for x in range(size[0]):
        for y in range(size[1]):
            for z in range(size[2]):
                val = scalars[x][y][z]
                val = val * val
                totalSum = totalSum + val
                closest = segment_array[x][y][z]
                chargePerAtom[closest] = chargePerAtom[closest] + val
                atomVolume[closest] = atomVolume[closest] + 1
    return chargePerAtom, atomVolume


def compute_atomic_charges(transitions, num_threads=4, same_atomic_positions=True, save_segmention=False):
    if (len(transitions) < 0):
        return
    data = transitions[0].hole_data
    segment_arrays = []
    if same_atomic_positions:
        segment_array = voronoi_segmentation(
            data.basis, data.atoms, data.size, num_threads)
        if save_segmention:
            segment_arrays.append(segment_array)
    for transition in transitions:
        data = transition.hole_data
        if not same_atomic_positions:
            segment_array = voronoi_segmentation(
                data.basis, data.atoms, data.size, num_threads)
            if save_segmention:
                segment_arrays.append(segment_array)
        transition.hole_charges, atomVolume = accumulate_atomic_charges(
            segment_array, data.scalars, data.atoms, data.size)
        data = transition.particle_data
        if not same_atomic_positions:
            segment_array = voronoi_segmentation(
                data.basis, data.atoms, data.size, num_threads)
            if save_segmention:
                segment_arrays.append(segment_array)
        transition.particle_charges, atomVolume = accumulate_atomic_charges(
            segment_array, data.scalars, data.atoms, data.size)
        transition.normalize()
    return segment_arrays


def read_atomic_charges(file_name):
    with open(file_name) as chargeFile:
        lines = chargeFile.readlines()
        hole_charges = []
        particle_charges = []
        total_hole_charge = 0
        total_particle_charge = 0
        for line in lines[1:]:
            splitted = line.split(",")
            hole_charge = float(splitted[3].strip())
            hole_charges.append(hole_charge)
            totalGSCharge = totalGSCharge + hole_charge
            particle_charge = float(splitted[4].strip())
            particle_charges.append(particle_charge)
            totalESCharge = totalESCharge + particle_charge
        for i in range(len(hole_charges)):
            hole_charges[i] = hole_charges[i] / total_hole_charge
        for i in range(len(particle_charges)):
            particle_charges[i] = particle_charges[i] / total_particle_charge
        return (hole_charges, particle_charges)


def save_atomic_charges(output_file, atoms, hole_charges, particle_charges, subgroup_names, atom_subgroup_map):
    with open(output_file, 'w') as f:
        f.write(
            "ID, Atomic number, Subgroup, Hole charge, Particle charge, Charge difference\n")
        for i in range(len(atoms)):
            f.write("%d, %d, %s, %f, %f, %f\n" % (
                atoms[i].id, atoms[i].atomic_number, subgroup_names[atom_subgroup_map[i]],
                hole_charges[i], particle_charges[i], particle_charges[i] - hole_charges[i]))


def read_metadata_file(csv_file: str):
    fileNames = []
    with open(csv_file, mode='r') as f:
        lines = f.readlines()
        for i in range(1, len(lines)):
            splitted = lines[i].strip().split(",")
            fileNames.append((splitted[0], splitted[1], splitted[2]))
    return fileNames


def load_subgroups(file_name: str):
    with open(file_name) as f:
        lines = f.readlines()
        splitted = lines[0].strip().split(",")
        subgroup_names = []
        for name in splitted:
            subgroup_names.append(name.strip())
        splitted = lines[1].strip().split(",")
        atom_subgroup_map = []
        for index in splitted:
            atom_subgroup_map.append(int(index.strip()))
        return (subgroup_names, atom_subgroup_map)


class SubgroupInfo:
    def __init__(self):
        self.init(0)

    def init(self, num_subgroups):
        self.num_subgroups = num_subgroups
        self.matrix = []
        self.subgroup_names = []
        self.hole_charges = []
        self.particle_charges = []
        self.atom_subgroup_map = []

    def set_subgroups(self, subgroup_names, atom_subgroup_map):
        self.init(len(subgroup_names))
        self.subgroup_names = subgroup_names
        self.atom_subgroup_map = atom_subgroup_map

    def __repr__(self):
        return ("Num groups: %d\nSubgroups: %s\nHole charges: %s\nParticle charges: %s\nMatrix:%s" %
                (self.num_subgroups, self.subgroup_names, self.hole_charges, self.particle_charges, self.matrix))

    def load_from_file(self, file_name: str):
        with open(file_name) as subgroupFile:
            lines = subgroupFile.readlines()
            names = lines[0].strip().split()
            self.init(len(names))
            for name in names:
                self.subgroup_names.append(name.upper())

            self.hole_charges = list(map(float, lines[1].strip().split()))
            self.particle_charges = list(map(float, lines[2].strip().split()))
            for i in range(self.num_subgroups):
                vals = list(map(float, lines[3 + i].strip().split()))
                self.matrix.append(vals)

    def save_to_file(self, file_name: str):
        with open(file_name, 'w') as f:
            for name in self.subgroup_names:
                f.write(name + " ")
            f.write("\n")
            for charge in self.hole_charges:
                f.write(str(charge) + " ")
            f.write("\n")
            for charge in self.particle_charges:
                f.write(str(charge) + " ")
            f.write("\n")
            for i in range(self.num_subgroups):
                for j in range(self.num_subgroups):
                    f.write("%f " % (self.matrix[i][j]))
                f.write("\n")


def _accumulate_subgroup_charges(hole_charges, particle_charges, num_subgroups, atom_subgroup_map):
    sumGSCharges = sum(hole_charges)
    sumESCharges = sum(particle_charges)
    ligandGSCharges = [0.0] * num_subgroups
    ligandESCharges = [0.0] * num_subgroups
    for i in range(len(atom_subgroup_map)):
        subgroupIndex = atom_subgroup_map[i]
        ligandGSCharges[subgroupIndex] = ligandGSCharges[subgroupIndex] + \
            hole_charges[i] / sumGSCharges
        ligandESCharges[subgroupIndex] = ligandESCharges[subgroupIndex] + \
            particle_charges[i] / sumESCharges
    return (ligandGSCharges, ligandESCharges)


def compute_subgroup_charges(transition, subgroup_info, use_hueristic=True):
    ligandGSCharges, ligandESCharges = _accumulate_subgroup_charges(transition.hole_charges, transition.particle_charges,
                                                                    subgroup_info.num_subgroups, subgroup_info.atom_subgroup_map)
    subgroup_info.hole_charges = ligandGSCharges
    subgroup_info.particle_charges = ligandESCharges
    if use_hueristic:
        chargeTransfer = _distribute_charges_heuristic(
            ligandGSCharges, ligandESCharges)
        subgroup_info.matrix = chargeTransfer
    else:
        chargeTransfer = _distribute_charges_optimize_quadratic(
            ligandGSCharges, ligandESCharges)
        subgroup_info.matrix = chargeTransfer


def _distribute_charges_heuristic(hole_charges, particle_charges):
    T_complete = []
    for i in range(len(hole_charges)):
        T_complete.append([0] * len(hole_charges))

    donors = []
    acceptors = []
    chargeDiff = []
    for i in range(len(hole_charges)):
        gsCharge = hole_charges[i]
        esCharge = particle_charges[i]
        if gsCharge > esCharge:
            donors.append(i)
        else:
            acceptors.append(i)
        diff = esCharge - gsCharge
        T_complete[i][i] = min(gsCharge, esCharge)
        chargeDiff.append(diff)

    totalAcceptorCharge = 0
    for acceptor in acceptors:
        totalAcceptorCharge = totalAcceptorCharge + chargeDiff[acceptor]

    # heuristic
    for donor in donors:
        chargeDeficit = -chargeDiff[donor]
        for acceptor in acceptors:
            contrib = chargeDeficit * \
                (chargeDiff[acceptor]) / totalAcceptorCharge
            T_complete[acceptor][donor] = contrib

    return T_complete


def _distribute_charges_optimize_quadratic(hole_charges, particle_charges):
    T_complete = []
    for i in range(len(hole_charges)):
        T_complete.append([0] * len(hole_charges))

    donors = []
    acceptors = []
    chargeDiff = []
    for i in range(len(hole_charges)):
        gsCharge = hole_charges[i]
        esCharge = particle_charges[i]
        if gsCharge > esCharge:
            donors.append(i)
        else:
            acceptors.append(i)
        diff = esCharge - gsCharge
        T_complete[i][i] = min(gsCharge, esCharge)
        chargeDiff.append(diff)

    totalAcceptorCharge = 0
    for acceptor in acceptors:
        totalAcceptorCharge = totalAcceptorCharge + chargeDiff[acceptor]

    n = len(donors)
    m = len(acceptors)
    if n == 0 or m == 0:
        return T_complete

    C = []
    d = []
    b = [totalAcceptorCharge / (n * m)] * (n * m)
    index = 0
    for donor in donors:
        chargeDeficit = -chargeDiff[donor]
        row = [0] * (n * m)
        for i in range(m):
            row[m * index + i] = 1
        C.append(row)
        d.append(chargeDeficit)
        index = index + 1

    index = 0
    for acceptor_index in range(m - 1):
        acceptor = acceptors[acceptor_index]
        chargeDeficit = chargeDiff[acceptor]
        row = [0] * (n * m)
        for i in range(n):
            row[m * i + index] = 1
        C.append(row)
        d.append(chargeDeficit)
        index = index + 1

    import numpy

    d = numpy.array(d)
    b = numpy.array(b)
    C = numpy.array(C)

    import cvxpy as cp

    # Construct the problem.
    x = cp.Variable(n * m)
    objective = cp.Minimize(cp.sum_squares(x - b))
    constraints = [x >= 0, C @ x == d]
    prob = cp.Problem(objective, constraints)

    # The optimal objective value is returned by `prob.solve()`.
    result = prob.solve()
    # The optimal value for x is stored in `x.value`.
    x_optimal = x.value

    index = 0
    for donor in donors:
        for acceptor in acceptors:
            contrib = x_optimal[index]
            T_complete[acceptor][donor] = contrib
            index = index + 1

    return T_complete


def create_diagram(subgroup_info, title="", show_plot=True, save_plot=False, file_name="diagram.pdf"):

    WIDTH = 700
    HEIGHT = 700
    PADDING_TOP = 0.05
    PADDING_LEFT = 0.08
    BAR_THICKNESS = 0.04
    GAP = 0.035
    FLOW_GAP = 0.0
    TEXT_GAP = 0.06
    TITLE_GAP = 0.12

    import matplotlib.pyplot as plt
    import matplotlib.path as mpath
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots()

    topPadding = HEIGHT * PADDING_TOP
    leftPadding = WIDTH * PADDING_LEFT
    barHeight = HEIGHT * BAR_THICKNESS
    textHeight = HEIGHT * TEXT_GAP
    titleHeight = HEIGHT * TITLE_GAP

    availableWidth = (1.0 - 2 * PADDING_LEFT -
                      (subgroup_info.num_subgroups - 1) * GAP) * WIDTH

    colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
              "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"]

    # draw ES bars
    xCoordsES = [0] * subgroup_info.num_subgroups
    total = 0.0
    topY = HEIGHT - topPadding - textHeight - titleHeight - barHeight

    plt.text(WIDTH / 2, HEIGHT - titleHeight, title,
             ha="center", family="serif", size=16)

    for i in range(subgroup_info.num_subgroups):
        x = leftPadding + total * availableWidth + i * GAP * WIDTH
        barWidth = subgroup_info.particle_charges[i] * availableWidth

        percent = "%.2f%%" % (subgroup_info.particle_charges[i] * 100)
        plt.text(x + barWidth / 2, topY + 6 + barHeight,
                 subgroup_info.subgroup_names[i], ha="center", family="sans-serif", size=10)
        plt.text(x + barWidth / 2, topY + 28 + barHeight, percent,
                 ha="center", family="monospace", size=10)
        rect = mpatches.Rectangle(
            [x, topY], barWidth, barHeight, faceColor=colors[i % len(colors)], ec="black", lw=0.5)
        ax.add_patch(rect)

        xCoordsES[i] = x
        total = total + subgroup_info.particle_charges[i]

    # draw GS bars
    xCoordsGS = [0] * subgroup_info.num_subgroups
    total = 0.0
    bottomY = topPadding + textHeight

    for i in range(subgroup_info.num_subgroups):
        x = leftPadding + total * availableWidth + i * GAP * WIDTH
        barWidth = subgroup_info.hole_charges[i] * availableWidth

        percent = "%.2f%%" % (subgroup_info.hole_charges[i] * 100)
        plt.text(x + barWidth / 2, bottomY - 22,
                 subgroup_info.subgroup_names[i], ha="center", family="sans-serif", size=10)
        plt.text(x + barWidth / 2, bottomY - 44, percent,
                 ha="center", family="monospace", size=10)
        rect = mpatches.Rectangle(
            [x, bottomY], barWidth, barHeight, faceColor=colors[i % len(colors)], ec="black", lw=0.5)
        ax.add_patch(rect)

        xCoordsGS[i] = x
        total = total + subgroup_info.hole_charges[i]

    # draw connectors
    propoFilledES = [0] * subgroup_info.num_subgroups
    propoFilledGS = [0] * subgroup_info.num_subgroups

    for i in range(subgroup_info.num_subgroups):
        xGS = xCoordsGS[i]
        for j in range(subgroup_info.num_subgroups):
            if subgroup_info.matrix[j][i] == 0.0:
                continue

            xES = xCoordsES[j]
            flow = subgroup_info.matrix[j][i]
            propGS = propoFilledGS[i]
            propES = propoFilledES[j]

            bottomXleft = xGS + availableWidth * propGS
            bottomXright = xGS + availableWidth * (propGS + flow)
            yb = bottomY + barHeight + FLOW_GAP * HEIGHT
            topXleft = xES + availableWidth * propES
            topXright = xES + availableWidth * (propES + flow)
            yt = topY - FLOW_GAP * HEIGHT
            yMiddle = (yt + yb) / 2

            # add a path patch
            Path = mpath.Path
            path_data = [
                (Path.MOVETO, [bottomXleft, yb]),
                (Path.CURVE4, [bottomXleft, yMiddle]),
                (Path.CURVE4, [topXleft, yMiddle]),
                (Path.CURVE4, [topXleft, yt]),
                (Path.LINETO, [topXright, yt]),
                (Path.CURVE4, [topXright, yMiddle]),
                (Path.CURVE4, [bottomXright, yMiddle]),
                (Path.CURVE4, [bottomXright, yb]),
                (Path.CLOSEPOLY, [bottomXleft, yb])
            ]
            codes, verts = zip(*path_data)
            path = mpath.Path(verts, codes)
            patch = mpatches.PathPatch(
                path, faceColor=colors[i % len(colors)], alpha=0.4, ec="none")
            ax.add_patch(patch)

            propoFilledGS[i] = propoFilledGS[i] + flow
            propoFilledES[j] = propoFilledES[j] + flow

    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(file_name, bbox_inches='tight')
