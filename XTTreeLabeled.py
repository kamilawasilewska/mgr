#
#    <CustomTools>
#      <Menu>
#       <Item name="XTTreeLabeled" icon="Python">
#         <Command>PythonXT::XT()</Command>
#       </Item>
#      </Menu>
#    </CustomTools>

import ImarisLib
import numpy as np


def XTTreeLabeled(imarisId):
    # Create an ImarisLib object
    imarisLib = ImarisLib.ImarisLib()

    # Get an imaris object with id aImarisId
    imarisApp = imarisLib.GetApplication(imarisId)

    # Get the currently loaded dataset
    currentImage = imarisApp.GetImage(0)

    # Get Cell data from Imaris
    cells = imarisApp.GetFactory().ToSpots(imarisApp.GetSurpassSelection())  # Works only for Spots mode

    cellsIds = cells.GetIds()  # ID cell form Statistic->Selection
    cellsTrackIds = cells.GetTrackIds()  # TrackID cell form Statistic->Selection
    cellTrackEdges = cells.GetTrackEdges()

    # Connection between TrackID and edges. Special definitions in function GetTrackIds.
    edges = np.column_stack((cellTrackEdges, cellsTrackIds))
    dataExtracted = extractUniqueData(edges)
    dataExtracted = np.column_stack((cellsIds, dataExtracted))  # Define array with CellID, TrackID, Egdes.

    # Get statistics from Imaris
    statistics = cells.GetStatistics()
    # uniqueSpecificName = np.transpose(np.unique(np.array(statistics.mNames)))
    data = np.column_stack((statistics.mIds, statistics.mNames, statistics.mValues))
    generations = data[np.where(data[:, 1] == 'Generation')]
    timeIndex = data[np.where(data[:, 1] == 'Time Index')]

    # Create table with CellId, TrackID, Generations, Time point
    table = prepareTableWithData(dataExtracted, generations, timeIndex)

    createUniqueLabels(table)
    res = np.array(table)
    a = 1

    new_names = ['Tree Labeled'] * len(table)  # aName
    new_values = [element[2] for element in table]  # aValue
    new_units = [''] * len(table)  # aUnits
    new_ids = cellsIds
    new_factor_names = ['Category', 'TreeLabeled']  # aFactorNames
    new_factors = [['Spot'] * len(table), [element[4] for element in table]]  # aFactor
    currentElement = imarisApp.GetSurpassSelection()

    currentElement.RemoveStatistics(new_names)
    currentElement.AddStatistics(new_names, new_values, new_units, new_factors, new_factor_names, new_ids)


def extractUniqueData(edges):
    dataExtracted = list()
    for cell in edges:
        if not [cell[2], cell[0]] in dataExtracted:
            dataExtracted.append([cell[2], cell[0]])
        if not [cell[2], cell[1]] in dataExtracted:
            dataExtracted.append([cell[2], cell[1]])
    dataExtracted.sort(key=lambda x: x[1])
    dataExtracted = np.array(dataExtracted)
    return dataExtracted


def prepareTableWithData(dataExtracted, generations, timeIndex):
    """
   :return [...[CellId, TrackID, Generations, TimePoint, LabelPlaceholder, CellGenerationIndex]]
    """
    table = []
    for i in range(len(dataExtracted)):
        if dataExtracted[i, 0] == int(generations[i, 0]):
            table.append(
                [
                    dataExtracted[i, 0],
                    dataExtracted[i, 1],
                    int(float(generations[i, 2])),
                    int(float(timeIndex[i, 2])),
                    None,
                    None
                ]
            )
    return table


def createUniqueLabels(table):
    """
   labels unique cells and adds number for cell during division
   :return void
    """
    for timeSlice in range(int(max(table, key=lambda x: x[3])[3])):
        imarisTimeSlice = timeSlice + 1
        currentSliceCells = [l for l in table if int(l[3]) == imarisTimeSlice]
        currentSliceCells = sorted(currentSliceCells, key=lambda x: x[2])
        if imarisTimeSlice == 1:
            initializeLabels(currentSliceCells)
            previousIteration = currentSliceCells
            continue
        for cell in currentSliceCells:
            thisIterationSameLabelCells = [l for l in currentSliceCells if int(l[1]) == cell[1]]
            prevIterationSameLabelCells = [l for l in previousIteration if int(l[1]) == cell[1]]
            if len(prevIterationSameLabelCells) >= len(thisIterationSameLabelCells):
                rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, table)
            elif len(prevIterationSameLabelCells) < len(thisIterationSameLabelCells):
                rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, table)
                try:
                    labeledCells
                except NameError:
                    labeledCells = 0
                try:
                    usedLabels
                except NameError:
                    usedLabels = []
                try:
                    label
                except NameError:
                    label = 1
                label, labeledCells, usedLabels = createNewLabel(
                    cell,
                    imarisTimeSlice,
                    prevIterationSameLabelCells,
                    previousIteration,
                    table,
                    thisIterationSameLabelCells,
                    label,
                    labeledCells,
                    usedLabels
                )
        try:
            del label, labeledCells, usedLabels
        except NameError:
            pass
        previousIteration = currentSliceCells


def createNewLabel(
        cell,
        imarisTimeSlice,
        prevIterationSameLabelCells,
        previousIteration,
        table,
        thisIterationSameLabelCells,
        label,
        labeledCells,
        usedLabels
):
    if not cell[1] in [cell[1] for cell in previousIteration]:
        labeledCells = findLabeledCellInGeneralTable(
            cell,
            imarisTimeSlice,
            labeledCells,
            table,
            thisIterationSameLabelCells
        )
    else:
        if label is None or label > 2:
            label = 1
        prevIterationCell = getParentCell(cell, prevIterationSameLabelCells, thisIterationSameLabelCells)
        createLabel(cell, prevIterationCell, label, usedLabels)
        cell[5] = labeledCells
        labeledCells += 1
        label += 1

    return label, labeledCells, usedLabels


def getParentCell(cell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    if len(prevIterationSameLabelCells) == 1:
        prevIterationCell = prevIterationSameLabelCells[0]
    else:
        prevIterationCell = [item for item in prevIterationSameLabelCells if
                             [item[1], item[2]] == [cell[1], cell[2]]]
    if not prevIterationCell:
        prevIterationCell = findParentCell(
            prevIterationCell,
            prevIterationSameLabelCells,
            thisIterationSameLabelCells
        )
    elif len(prevIterationCell) > 1:
        for prev in prevIterationCell:
            if not isinstance(prev, list):
                break
            if len(prev[4].split('_')) - 1 == cell[2]:
                prevIterationCell = prev
                break
    return prevIterationCell


def findLabeledCellInGeneralTable(cell, imarisTimeSlice, labeledCells, table, thisIterationSameLabelCells):
    tmp = allPreviouslyLabeledCells(imarisTimeSlice, table)
    for labeledCell in tmp:
        if labeledCell[1] == cell[1] and not labeledCell[5] in [item[-1] for item in
                                                                thisIterationSameLabelCells]:
            cell[4] = labeledCell[4]
            cell[5] = labeledCells
            labeledCells += 1
            break
    del tmp
    return labeledCells


def findParentCell(prevIterationCell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    for thisIterationCell in thisIterationSameLabelCells:
        for prev in prevIterationSameLabelCells:
            if prev[2] != thisIterationCell[2] and prev[4] not in [e[4] for e in thisIterationSameLabelCells]:
                prevIterationCell = prev
                break
        if prevIterationCell:
            break
    if not prevIterationCell:
        for parentCell in prevIterationSameLabelCells:
            if len([e[2] for e in thisIterationSameLabelCells if e[4] == parentCell[4]]) > 1:
                continue

            prevIterationCell = parentCell
            break
    return prevIterationCell


def createLabel(cell, prevIterationCell, label, usedLabels):
    if any(isinstance(i, list) for i in prevIterationCell):
        prevIterationCell = prevIterationCell[0]
    newLabel = prevIterationCell[4]
    if len(str(newLabel).split('_')) - 1 < cell[2]:
        newLabel += '_'
        newLabel += str(label)
        label += 1
        usedLabels.append(prevIterationCell[4])
    cell[4] = newLabel
    return label


def initializeLabels(currentSliceCells):
    label = 1
    for cell in currentSliceCells:
        cell[4] = str(label)
        cell[5] = True
        label += 1
    return label


def allPreviouslyLabeledCells(imarisTimeSlice, table):
    tmp = []
    slice_length = len([l for l in table if l[3] < imarisTimeSlice]) - 1
    for i in range(slice_length):
        tmp.append(table[slice_length - i])
    return tmp


def rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, table):
    if len(prevIterationSameLabelCells) != len(thisIterationSameLabelCells) and \
            [cell[0], cell[1], cell[2], cell[3]] not in [
        [
            element[0],
            element[1],
            element[2],
            element[3]
        ] for element in allPreviouslyLabeledCells(cell[3], table)]:
        return
    for labeledCell in prevIterationSameLabelCells:
        if not [labeledCell[1], labeledCell[5]] in [[item[1], item[-1]] for item in thisIterationSameLabelCells]:
            cell[4] = labeledCell[4]
            cell[5] = labeledCell[5]
            break
