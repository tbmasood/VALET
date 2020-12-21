def write_segmentation(output_file, segment_array, basis, array_name="density"):
    import vtk
    import numpy
    from vtk.numpy_interface import dataset_adapter as dsa
    
    size = segment_array.shape
    segment_array = segment_array.flatten('F')

    grid = vtk.vtkImageData()
    grid.SetOrigin(0, 0, 0)
    grid.SetSpacing(basis[0][0], basis[1][1], basis[2][2])
    grid.SetDimensions(size)
    grid.GetPointData().SetScalars(dsa.numpyTovtkDataArray(segment_array, array_name))

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(grid)
    writer.Write()


def write_atoms(output_file, atoms, atomic_charges):
    import vtk
    import numpy
    import electrans.atomic_data as atomic_data

    pointsArr = vtk.vtkPoints()
    atomTypesArr = vtk.vtkIntArray()
    atomTypesArr.SetName("atomType")
    atomRadiusArr = vtk.vtkFloatArray()
    atomRadiusArr.SetName("atomRadius")
    atomColorArr = vtk.vtkFloatArray()
    atomColorArr.SetName("atomColor")
    atomColorArr.SetNumberOfComponents(3)
    atomChargeArr = vtk.vtkFloatArray()
    atomChargeArr.SetName("charge")
    for atom in atoms:
        x = atom.coordinates[0]
        y = atom.coordinates[1]
        z = atom.coordinates[2]
        pointsArr.InsertNextPoint(x, y, z)
        atomTypesArr.InsertNextTuple1(atom.atomic_number)
        radius = atom.radius
        atomRadiusArr.InsertNextTuple1(radius)
        color = atomic_data.color(atomic_data.atomic_symbol(atom.atomic_number))
        atomColorArr.InsertNextTuple3(color[0], color[1], color[2])
        atomChargeArr.InsertNextTuple1(atomic_charges[atom.id])

    # Add bonds
    bondsArr = vtk.vtkCellArray()
    for i in range(len(atoms)):
        xi = atoms[i].coordinates[0]
        yi = atoms[i].coordinates[1]
        zi = atoms[i].coordinates[2]
        radius_i = atoms[i].radius
        for j in range(i+1, len(atoms)):
            xj = atoms[j].coordinates[0]
            yj = atoms[j].coordinates[1]
            zj = atoms[j].coordinates[2]
            radius_j = atoms[j].radius
            xDiff = xi - xj
            yDiff = yi - yj
            zDiff = zi - zj
            dist = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff
            if dist < 0.4 * (radius_i + radius_j) * (radius_i + radius_j):
                bond = vtk.vtkLine()
                bond.GetPointIds().SetId(0, i)
                bond.GetPointIds().SetId(1, j)
                bondsArr.InsertNextCell(bond)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(pointsArr)
    polydata.SetLines(bondsArr)
    polydata.GetPointData().AddArray(atomTypesArr)
    polydata.GetPointData().AddArray(atomRadiusArr)
    polydata.GetPointData().AddArray(atomColorArr)
    polydata.GetPointData().AddArray(atomChargeArr)
    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(polydata)
    writer.Write()


def write_segments(output_file, segment_array, basis, atoms):
    import vtk
    import numpy
    from vtk.numpy_interface import dataset_adapter as dsa
    from scipy.ndimage.filters import gaussian_filter
    
    size = segment_array.shape
    mask = numpy.zeros(size)
    for x in range(size[0]):
        for y in range(size[1]):
            for z in range(size[2]):
                pt = (x * basis[0][0], y * basis[1][1], z * basis[2][2])
                for i in range(len(atoms)):
                    radius = atoms[i].radius
                    xDiff = atoms[i].coordinates[0] - pt[0]
                    yDiff = atoms[i].coordinates[1] - pt[1]
                    zDiff = atoms[i].coordinates[2] - pt[2]
                    dist = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff
                    if dist < 4 * radius * radius:
                        mask[x][y][z] = 1
                        break

    grid = vtk.vtkImageData()
    grid.SetOrigin(0, 0, 0)
    grid.SetSpacing(basis[0][0], basis[1][1], basis[2][2])
    grid.SetDimensions(size)

    appendedData = vtk.vtkAppendPolyData()

    for i in range(len(atoms)):
        selected = numpy.where(segment_array == i, mask, 0)
        selected = gaussian_filter(selected, sigma=0.4)
        selected = selected.flatten('F')
        grid.GetPointData().SetScalars(dsa.numpyTovtkDataArray(selected, "seg"))

        surface = vtk.vtkMarchingCubes()
        surface.SetInputData(grid)
        surface.ComputeNormalsOn()
        surface.SetValue(0, 0.4)
        surface.Update()
        scalarArray = surface.GetOutput().GetPointData().GetScalars()
        for j in range(surface.GetOutput().GetNumberOfPoints()):
            scalarArray.SetTuple1(j, i)

        appendedData.AddInputData(surface.GetOutput())
        appendedData.Update()

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendedData.GetOutputPort())
    cleanFilter.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(cleanFilter.GetOutput())
    writer.Write()


def write_subgroup_segments(output_file, segment_array, basis, atoms, ligands):
    import vtk
    import numpy
    from vtk.numpy_interface import dataset_adapter as dsa
    from scipy.ndimage.filters import gaussian_filter
    
    size = segment_array.shape
    mask = numpy.zeros(size)
    for x in range(size[0]):
        for y in range(size[1]):
            for z in range(size[2]):
                pt = (x * basis[0][0], y * basis[1][1], z * basis[2][2])
                for i in range(len(atoms)):
                    radius = atoms[i].radius
                    xDiff = atoms[i].coordinates[0] - pt[0]
                    yDiff = atoms[i].coordinates[1] - pt[1]
                    zDiff = atoms[i].coordinates[2] - pt[2]
                    dist = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff
                    if dist < 4 * radius * radius:
                        mask[x][y][z] = 1
                        break

    grid = vtk.vtkImageData()
    grid.SetOrigin(0, 0, 0)
    grid.SetSpacing(basis[0][0], basis[1][1], basis[2][2])
    grid.SetDimensions(size)

    appendedData = vtk.vtkAppendPolyData()

    for i in range(len(ligands)):
        selected = numpy.zeros(size)
        for atomID in ligands[i]:
            selected = selected + numpy.where(segment_array == atomID, mask, 0)

        selected = gaussian_filter(selected, sigma=0.4)
        selected = selected.flatten('F')
        grid.GetPointData().SetScalars(dsa.numpyTovtkDataArray(selected, "seg"))

        surface = vtk.vtkMarchingCubes()
        surface.SetInputData(grid)
        surface.ComputeNormalsOn()
        surface.SetValue(0, 0.4)
        surface.Update()
        scalarArray = surface.GetOutput().GetPointData().GetScalars()
        for j in range(surface.GetOutput().GetNumberOfPoints()):
            scalarArray.SetTuple1(j, i)

        appendedData.AddInputData(surface.GetOutput())
        appendedData.Update()

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendedData.GetOutputPort())
    cleanFilter.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(cleanFilter.GetOutput())
    writer.Write()
