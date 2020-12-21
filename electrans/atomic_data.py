_elements = [
    "Xx", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og"
]

_default_atom = {'rgb': [0, 0, 0], 'vdwRadius': 3.00, 'covalentRadius': 3.00}

_atom_colors = {
    "Xx": _default_atom,
    "H": {'rgb': [255, 255, 255], 'vdwRadius': 1.10, 'covalentRadius': 0.31},
    "He": {'rgb': [217, 255, 255], 'vdwRadius': 1.40, 'covalentRadius': 0.28},
    "Li": {'rgb': [204, 128, 255], 'vdwRadius': 1.81, 'covalentRadius': 1.28},
    "Be": {'rgb': [194, 255, 0], 'vdwRadius': 1.53, 'covalentRadius': 0.96},
    "B": {'rgb': [255, 181, 181], 'vdwRadius': 1.92, 'covalentRadius': 0.84},
    "C": {'rgb': [144, 144, 144], 'vdwRadius': 1.70, 'covalentRadius': 0.76},
    "N": {'rgb': [48, 80, 248], 'vdwRadius': 1.55, 'covalentRadius': 0.71},
    "O": {'rgb': [255, 13, 13], 'vdwRadius': 1.52, 'covalentRadius': 0.66},
    "F": {'rgb': [144, 224, 80], 'vdwRadius': 1.47, 'covalentRadius': 0.57},
    "Ne": {'rgb': [179, 227, 245], 'vdwRadius': 1.54, 'covalentRadius': 0.58},
    "Na": {'rgb': [171, 92, 242], 'vdwRadius': 2.27, 'covalentRadius': 1.66},
    "Mg": {'rgb': [138, 255, 0], 'vdwRadius': 1.73, 'covalentRadius': 1.41},
    "Al": {'rgb': [191, 166, 166], 'vdwRadius': 1.84, 'covalentRadius': 1.21},
    "Si": {'rgb': [240, 200, 160], 'vdwRadius': 2.10, 'covalentRadius': 1.11},
    "P": {'rgb': [255, 128, 0], 'vdwRadius': 1.80, 'covalentRadius': 1.07},
    "S": {'rgb': [255, 255, 48], 'vdwRadius': 1.80, 'covalentRadius': 1.05},
    "Cl": {'rgb': [31, 240, 31], 'vdwRadius': 1.75, 'covalentRadius': 1.02},
    "Ar": {'rgb': [128, 209, 227], 'vdwRadius': 1.88, 'covalentRadius': 1.06},
    "K": {'rgb': [143, 64, 212], 'vdwRadius': 2.75, 'covalentRadius': 2.03},
    "Ca": {'rgb': [61, 255, 0], 'vdwRadius': 2.31, 'covalentRadius': 1.76},
    "Sc": {'rgb': [230, 230, 230], 'vdwRadius': 2.30, 'covalentRadius': 1.70},
    "Ti": {'rgb': [191, 194, 199], 'vdwRadius': 2.15, 'covalentRadius': 1.60},
    "V": {'rgb': [166, 166, 171], 'vdwRadius': 2.05, 'covalentRadius': 1.53},
    "Cr": {'rgb': [138, 153, 199], 'vdwRadius': 2.05, 'covalentRadius': 1.39},
    "Mn": {'rgb': [156, 122, 199], 'vdwRadius': 2.05, 'covalentRadius': 1.39},
    "Fe": {'rgb': [224, 102, 51], 'vdwRadius': 2.05, 'covalentRadius': 1.32},
    "Co": {'rgb': [240, 144, 160], 'vdwRadius': 2.00, 'covalentRadius': 1.26},
    "Ni": {'rgb': [80, 208, 80], 'vdwRadius': 2.00, 'covalentRadius': 1.24},
    "Cu": {'rgb': [200, 128, 51], 'vdwRadius': 2.00, 'covalentRadius': 1.32},
    "Zn": {'rgb': [125, 128, 176], 'vdwRadius': 2.10, 'covalentRadius': 1.22},
    "Ga": {'rgb': [194, 143, 143], 'vdwRadius': 1.87, 'covalentRadius': 1.22},
    "Ge": {'rgb': [102, 143, 143], 'vdwRadius': 2.11, 'covalentRadius': 1.20},
    "As": {'rgb': [189, 128, 227], 'vdwRadius': 1.85, 'covalentRadius': 1.19},
    "Se": {'rgb': [255, 161, 0], 'vdwRadius': 1.90, 'covalentRadius': 1.20},
    "Br": {'rgb': [166, 41, 41], 'vdwRadius': 1.83, 'covalentRadius': 1.20},
    "Kr": {'rgb': [92, 184, 209], 'vdwRadius': 2.02, 'covalentRadius': 1.16},
    "Rb": {'rgb': [112, 46, 176], 'vdwRadius': 3.03, 'covalentRadius': 1.20},
    "Sr": {'rgb': [0, 255, 0], 'vdwRadius': 2.49, 'covalentRadius': 1.95},
    "Y": {'rgb': [148, 255, 255], 'vdwRadius': 2.40, 'covalentRadius': 1.90},
    "Zr": {'rgb': [148, 224, 224], 'vdwRadius': 2.30, 'covalentRadius': 1.75},
    "Nb": {'rgb': [115, 194, 201], 'vdwRadius': 2.15, 'covalentRadius': 1.64},
    "Mo": {'rgb': [84, 181, 181], 'vdwRadius': 2.10, 'covalentRadius': 1.54},
    "Tc": {'rgb': [59, 158, 158], 'vdwRadius': 2.05, 'covalentRadius': 1.47},
    "Ru": {'rgb': [36, 143, 143], 'vdwRadius': 2.05, 'covalentRadius': 1.46},
    "Rh": {'rgb': [10, 125, 140], 'vdwRadius': 2.00, 'covalentRadius': 1.42},
    "Pd": {'rgb': [0, 105, 133], 'vdwRadius': 2.05, 'covalentRadius': 1.39},
    "Ag": {'rgb': [192, 192, 192], 'vdwRadius': 2.10, 'covalentRadius': 1.45},
    "Cd": {'rgb': [255, 217, 143], 'vdwRadius': 2.20, 'covalentRadius': 1.44},
    "In": {'rgb': [166, 117, 115], 'vdwRadius': 2.20, 'covalentRadius': 1.42},
    "Sn": {'rgb': [102, 128, 128], 'vdwRadius': 1.93, 'covalentRadius': 1.39},
    "Sb": {'rgb': [158, 99, 181], 'vdwRadius': 2.17, 'covalentRadius': 1.39},
    "Te": {'rgb': [212, 122, 0], 'vdwRadius': 2.06, 'covalentRadius': 1.38},
    "I": {'rgb': [148, 0, 148], 'vdwRadius': 1.98, 'covalentRadius': 1.39},
    "Xe": {'rgb': [66, 158, 176], 'vdwRadius': 2.16, 'covalentRadius': 1.40},
    "Cs": {'rgb': [87, 23, 143], 'vdwRadius': 3.43, 'covalentRadius': 2.44},
    "Ba": {'rgb': [0, 201, 0], 'vdwRadius': 2.68, 'covalentRadius': 2.15},
    "La": {'rgb': [112, 212, 255], 'vdwRadius': 2.50, 'covalentRadius': 2.07},
    "Ce": {'rgb': [255, 255, 199], 'vdwRadius': 2.48, 'covalentRadius': 2.04},
    "Pr": {'rgb': [217, 255, 199], 'vdwRadius': 2.47, 'covalentRadius': 2.03},
    "Nd": {'rgb': [199, 255, 199], 'vdwRadius': 2.45, 'covalentRadius': 2.01},
    "Pm": {'rgb': [163, 255, 199], 'vdwRadius': 2.43, 'covalentRadius': 1.99},
    "Sm": {'rgb': [143, 255, 199], 'vdwRadius': 2.42, 'covalentRadius': 1.98},
    "Eu": {'rgb': [97, 255, 199], 'vdwRadius': 2.40, 'covalentRadius': 1.98},
    "Gd": {'rgb': [69, 255, 199], 'vdwRadius': 2.38, 'covalentRadius': 1.96},
    "Tb": {'rgb': [48, 255, 199], 'vdwRadius': 2.37, 'covalentRadius': 1.94},
    "Dy": {'rgb': [31, 255, 199], 'vdwRadius': 2.35, 'covalentRadius': 1.92},
    "Ho": {'rgb': [0, 255, 156], 'vdwRadius': 2.33, 'covalentRadius': 1.92},
    "Er": {'rgb': [0, 230, 117], 'vdwRadius': 2.32, 'covalentRadius': 1.89},
    "Tm": {'rgb': [0, 212, 82], 'vdwRadius': 2.30, 'covalentRadius': 1.90},
    "Yb": {'rgb': [0, 191, 56], 'vdwRadius': 2.28, 'covalentRadius': 1.87},
    "Lu": {'rgb': [0, 171, 36], 'vdwRadius': 2.27, 'covalentRadius': 1.75},
    "Hf": {'rgb': [77, 194, 255], 'vdwRadius': 2.25, 'covalentRadius': 1.87},
    "Ta": {'rgb': [77, 166, 255], 'vdwRadius': 2.20, 'covalentRadius': 1.70},
    "W": {'rgb': [33, 148, 214], 'vdwRadius': 2.10, 'covalentRadius': 1.62},
    "Re": {'rgb': [38, 125, 171], 'vdwRadius': 2.05, 'covalentRadius': 1.51},
    "Os": {'rgb': [38, 102, 150], 'vdwRadius': 2.00, 'covalentRadius': 1.44},
    "Ir": {'rgb': [23, 84, 135], 'vdwRadius': 2.00, 'covalentRadius': 1.41},
    "Pt": {'rgb': [208, 208, 224], 'vdwRadius': 2.05, 'covalentRadius': 1.36},
    "Au": {'rgb': [255, 209, 35], 'vdwRadius': 2.10, 'covalentRadius': 1.36},
    "Hg": {'rgb': [184, 184, 208], 'vdwRadius': 2.05, 'covalentRadius': 1.32},
    "Tl": {'rgb': [166, 84, 77], 'vdwRadius': 1.96, 'covalentRadius': 1.45},
    "Pb": {'rgb': [87, 89, 97], 'vdwRadius': 2.02, 'covalentRadius': 1.46},
    "Bi": {'rgb': [158, 79, 181], 'vdwRadius': 2.07, 'covalentRadius': 1.48},
    "Po": {'rgb': [171, 92, 0], 'vdwRadius': 1.97, 'covalentRadius': 1.40},
    "At": {'rgb': [117, 79, 69], 'vdwRadius': 2.02, 'covalentRadius': 1.50},
    "Rn": {'rgb': [66, 130, 150], 'vdwRadius': 2.20, 'covalentRadius': 1.50},
    "Fr": {'rgb': [66, 0, 102], 'vdwRadius': 3.48, 'covalentRadius': 2.60},
    "Ra": {'rgb': [0, 125, 0], 'vdwRadius': 2.83, 'covalentRadius': 2.21},
    "Ac": {'rgb': [112, 171, 250], 'vdwRadius': 2.00, 'covalentRadius': 2.15},
    "Th": {'rgb': [0, 186, 255], 'vdwRadius': 2.40, 'covalentRadius': 2.06},
    "Pa": {'rgb': [0, 161, 255], 'vdwRadius': 2.00, 'covalentRadius': 2.00},
    "U": {'rgb': [0, 143, 255], 'vdwRadius': 2.30, 'covalentRadius': 1.96},
    "Np": {'rgb': [0, 128, 255], 'vdwRadius': 2.00, 'covalentRadius': 1.90},
    "Pu": {'rgb': [0, 107, 255], 'vdwRadius': 2.00, 'covalentRadius': 1.87},
    "Am": {'rgb': [84, 92, 242], 'vdwRadius': 2.00, 'covalentRadius': 1.80},
    "Cm": {'rgb': [120, 92, 227], 'vdwRadius': 2.00, 'covalentRadius': 1.69},
    "Bk": {'rgb': [138, 79, 227], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Cf": {'rgb': [161, 54, 212], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Es": {'rgb': [179, 31, 212], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Fm": {'rgb': [179, 31, 186], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Md": {'rgb': [179, 13, 166], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "No": {'rgb': [189, 13, 135], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Lr": {'rgb': [199, 0, 102], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Rf": {'rgb': [204, 0, 89], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Db": {'rgb': [209, 0, 79], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Sg": {'rgb': [217, 0, 69], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Bh": {'rgb': [224, 0, 56], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Hs": {'rgb': [230, 0, 46], 'vdwRadius': 2.00, 'covalentRadius': 1.70},
    "Mt": {'rgb': [235, 0, 38], 'vdwRadius': 2.00, 'covalentRadius': 1.70}
}


def atomic_symbol(atomic_number) -> str:
    if atomic_number >= len(_elements):
        return _elements[0]
    else:
        return _elements[atomic_number]


def color(name:str):
    rgb = _atom_colors.get(name, _default_atom)['rgb']
    rgb.append(255)
    return [rgb[0]/255.0, rgb[1]/255.0, rgb[2]/255.0, rgb[3]/255.0]


def vdw_radius(name:str):
    return _atom_colors.get(name, _default_atom)["vdwRadius"]


def covalent_radius(name:str):
    return _atom_colors.get(name, _default_atom)["covalentRadius"]
