class Atom:
    def __init__(self, id, coordinates, atomic_number):
        self.id = id
        self.coordinates = coordinates
        self.atomic_number = atomic_number
        self.radius = 0.0

    def __repr__(self):
        return ("%d (%d): %s" % (self.id, self.atomic_number, self.coordinates))


class CubefileData:
    def __init__(self):
        self.basis = []
        self.size = (0, 0, 0)
        self.scalars = []
        self.atoms = []

    def num_atoms(self):
        return len(self.atoms)


def load_cubefile(input_file: str):
    from pathlib import Path
    import numpy
    import electrans.atomic_data as atomic_data

    filePath = Path(input_file)
    if filePath.suffix == ".xz":
        import lzma
        with lzma.open(filePath, mode='rt') as f:
            lines = f.readlines()
    elif filePath.suffix == ".bzip2" or filePath.suffix == ".bz2":
        import bz2
        with bz2.open(filePath, mode='rt') as f:
            lines = f.readlines()
    elif filePath.suffix == ".gz":
        import gzip
        with gzip.open(filePath, mode='rt') as f:
            lines = f.readlines()
    else:
        with open(filePath, mode='r') as f:
            lines = f.readlines()

    cubefile_data = CubefileData()

    nVal = 1
    extraLines = 0
    splitted = lines[2].strip().split()
    numAtoms = int(splitted[0])
    origin = numpy.array(list(map(float, splitted[1:4])))
    if numAtoms < 0:
        numAtoms = -numAtoms
        extraLines = 1
    if len(splitted) > 4:
        nVal = int(splitted[4])

    bohrToAngstrom = 0.529177
    unitsInAngstrom = False
    splitted = lines[3].strip().split()
    sizeX = int(splitted[0])
    if sizeX < 0:
        sizeX = -sizeX
        unitsInAngstrom = True
    a1 = numpy.array(list(map(float, splitted[1:4])))

    splitted = lines[4].strip().split()
    sizeY = abs(int(splitted[0]))
    a2 = numpy.array(list(map(float, splitted[1:4])))

    splitted = lines[5].strip().split()
    sizeZ = abs(int(splitted[0]))
    a3 = numpy.array(list(map(float, splitted[1:4])))

    cubefile_data.size = (sizeX, sizeY, sizeZ)
    cubefile_data.basis = numpy.array([a1, a2, a3])
    if not unitsInAngstrom:
        cubefile_data.basis = bohrToAngstrom * cubefile_data.basis

    cubefile_data.atoms = []
    for i in range(numAtoms):
        splitted = lines[6 + i].strip().split()
        atomic_number = int(splitted[0])
        atomPos = numpy.array(list(map(float, splitted[2:5])))
        atomPos = atomPos - origin
        if not unitsInAngstrom:
            atomPos = bohrToAngstrom * atomPos
        atom = Atom(i, atomPos, atomic_number)
        atomRadius = atomic_data.vdw_radius(
            atomic_data.atomic_symbol(atomic_number))
        atom.radius = atomRadius
        cubefile_data.atoms.append(atom)

    chg = []
    for i in range(6 + numAtoms + extraLines, len(lines)):
        splitted = lines[i].strip().split()
        for val in splitted:
            chg.append(float(val))

    chgdata = numpy.array(chg).astype(numpy.float32)
    chosenIdx = 0
    chgdata = chgdata[chosenIdx::nVal]
    cubefile_data.scalars = chgdata.reshape(cubefile_data.size)

    return cubefile_data