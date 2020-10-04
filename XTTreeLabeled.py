#    <CustomTools>
#      <Menu>
#       <Item name="XTTreeLabeled" icon="Python">
#         <Command>PythonXT::XT()</Command>
#       </Item>
#      </Menu>
#    </CustomTools>

import ImarisLib
import numpy as np

CELL_LABEL_DELIMITER = '_'

CELL_ID = 0
CELL_LINE_NUMBER = 1
CELL_GENERATIONS = 2
CELL_TIME_SLICE = 3
CELL_LABEL = 4
CELLS_LABELED_TOTAL = 5
CELL_DIVIDED = 6


def XTTreeLabeled(imarisId):
    # Create an ImarisLib object
    imarisLib = ImarisLib.ImarisLib()

    # Get an imaris object with id aImarisId
    imarisApp = imarisLib.GetApplication(imarisId)

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
    uniqueSpecificName = np.transpose(np.unique(np.array(statistics.mNames)))
    data = np.column_stack((statistics.mIds, statistics.mNames, statistics.mValues))
    generations = data[np.where(data[:, 1] == 'Generation')]
    timeIndex = data[np.where(data[:, 1] == 'Time Index')]
    timeSincePreviousDivision = data[np.where(data[:, 1] == 'Time Since Previous Division')]

    # Create table with CellId, TrackID, Generations, Time point
    table = prepareTableWithData(dataExtracted, generations, timeIndex, timeSincePreviousDivision)

    createUniqueLabels(table)
    new_table = np.array(table)

    new_names = ['Tree Labeled'] * len(table)  # aName
    new_values = [element[CELL_GENERATIONS] for element in table]  # aValue
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
        if not [cell[CELL_GENERATIONS], cell[CELL_ID]] in dataExtracted:
            dataExtracted.append([cell[CELL_GENERATIONS], cell[CELL_ID]])
        if not [cell[CELL_GENERATIONS], cell[CELL_LINE_NUMBER]] in dataExtracted:
            dataExtracted.append([cell[CELL_GENERATIONS], cell[CELL_LINE_NUMBER]])
    dataExtracted.sort(key=lambda x: x[CELL_LINE_NUMBER])
    dataExtracted = np.array(dataExtracted)
    return dataExtracted


def prepareTableWithData(dataExtracted, generations, timeIndex, timeSincePreviousDivision):
    """
   :return [...[CELL_ID, CELL_LINE_NUMBER, CELL_GENERATIONS, CELL_TIME_SLICE, CELL_LABEL, CELLS_LABELED_TOTAL, CELL_DIVIDED]]
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
                    None,
                    None
                ]
            )

    cellsAfterDivision = []
    for i in range(len(timeSincePreviousDivision)):
        if timeSincePreviousDivision[i, 2] == '1.0':
            cellsAfterDivision.append([timeSincePreviousDivision[i, 0]])

    for j in range(len(table)):
        for i in range(len(cellsAfterDivision)):
            if int(cellsAfterDivision[i][CELL_ID]) == table[j][CELL_ID]:
                table[j][-1] = 1
                break
            else:
                table[j][-1] = 0
    return table


def createUniqueLabels(table):
    """
   labels unique cells and adds number for cell during division
   :return void
    """
    for timeSlice in range(int(max(table, key=lambda x: x[CELL_TIME_SLICE])[CELL_TIME_SLICE])):
        imarisTimeSlice = timeSlice + 1
        currentSliceCells = [l for l in table if int(l[CELL_TIME_SLICE]) == imarisTimeSlice]
        currentSliceCells = sorted(currentSliceCells, key=lambda x: x[CELL_GENERATIONS])
        if imarisTimeSlice == 1:
            initializeLabels(currentSliceCells)
            previousIteration = currentSliceCells
            continue
        for cell in currentSliceCells:
            thisIterationSameLabelCells = [l for l in currentSliceCells if int(l[CELL_LINE_NUMBER]) == cell[CELL_LINE_NUMBER]]
            prevIterationSameLabelCells = [l for l in previousIteration if int(l[CELL_LINE_NUMBER]) == cell[CELL_LINE_NUMBER]]

            rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, table)
            if not (cell[CELL_LABEL] is None or cell[CELLS_LABELED_TOTAL] is None):
                continue

            try:
                labeledCells
            except NameError:
                labeledCells = 0
            try:
                label
            except NameError:
                label = 1
            if cell[CELL_DIVIDED] == 0:
                labeledCells = rewriteFromAllLabeledCells(cell, imarisTimeSlice, labeledCells, table, thisIterationSameLabelCells)
            if not (cell[CELL_LABEL] is None or cell[CELLS_LABELED_TOTAL] is None):
                continue

            label, labeledCells = createNewLabel(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, label, labeledCells)
        try:
            del label, labeledCells
        except NameError:
            pass
        previousIteration = currentSliceCells


def createNewLabel(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, label, labeledCells):
    if label is None or label > 2:
        label = 1
    prevIterationCell = getParentCell(cell, prevIterationSameLabelCells, thisIterationSameLabelCells)
    createLabel(cell, prevIterationCell, label)
    cell[CELLS_LABELED_TOTAL] = labeledCells
    labeledCells += 1
    label += 1

    return label, labeledCells


def getParentCell(cell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    if len(prevIterationSameLabelCells) == CELL_LINE_NUMBER:
        prevIterationCell = prevIterationSameLabelCells[CELL_ID]
    else:
        prevIterationCell = []
        for item in prevIterationSameLabelCells:
            if item[CELL_LINE_NUMBER] == cell[CELL_LINE_NUMBER] and item[CELL_GENERATIONS] != cell[CELL_GENERATIONS]:
                prevIterationCell.append(item)

    if not prevIterationCell:
        prevIterationCell = findParentCell(prevIterationCell, prevIterationSameLabelCells, thisIterationSameLabelCells)
    elif len(prevIterationCell) > CELL_LINE_NUMBER:
        for prev in prevIterationCell:
            if not isinstance(prev, list):
                break
            if len(prev[CELL_LABEL].split(CELL_LABEL_DELIMITER)) - CELL_LINE_NUMBER == cell[CELL_GENERATIONS]:
                prevIterationCell = prev
                break
    return prevIterationCell


def rewriteFromAllLabeledCells(cell, imarisTimeSlice, labeledCells, table, thisIterationSameLabelCells):
    tmp = allPreviouslyLabeledCells(imarisTimeSlice, table)
    for labeledCell in tmp:
        if labeledCell[CELL_LINE_NUMBER] == cell[CELL_LINE_NUMBER] \
                and labeledCell[CELL_GENERATIONS] == cell[CELL_GENERATIONS] \
                and not labeledCell[CELL_LABEL] in [item[CELL_LABEL] for item in thisIterationSameLabelCells] \
                :
            cell[CELL_LABEL] = labeledCell[CELL_LABEL]
            cell[CELLS_LABELED_TOTAL] = labeledCells
            labeledCells += 1
            break
    return labeledCells


def findParentCell(prevIterationCell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    for thisIterationCell in thisIterationSameLabelCells:
        for prev in prevIterationSameLabelCells:
            if prev[CELL_GENERATIONS] != thisIterationCell[CELL_GENERATIONS] and prev[CELL_LABEL] not in [e[CELL_LABEL] for e in thisIterationSameLabelCells]:
                prevIterationCell = prev
                break
        if prevIterationCell:
            break
    if not prevIterationCell:
        for parentCell in prevIterationSameLabelCells:
            if len([e[CELL_GENERATIONS] for e in thisIterationSameLabelCells if e[CELL_LABEL] == parentCell[CELL_LABEL]]) > 1:
                continue

            prevIterationCell = parentCell
            break
    return prevIterationCell


def createLabel(cell, prevIterationCell, label):
    if any(isinstance(i, list) for i in prevIterationCell):
        prevIterationCell = prevIterationCell[CELL_ID]
    newLabel = prevIterationCell[CELL_LABEL]
    if len(str(newLabel).split('_')) - 1 < cell[CELL_GENERATIONS]:
        newLabel += '_'
        newLabel += str(label)
        label += 1
    cell[CELL_LABEL] = newLabel
    return label


def initializeLabels(currentSliceCells):
    label = 1
    for cell in currentSliceCells:
        cell[CELL_LABEL] = str(label)
        cell[CELLS_LABELED_TOTAL] = 1
        label += 1


def allPreviouslyLabeledCells(imarisTimeSlice, table):
    tmp = []
    slice_length = len([l for l in table if l[CELL_TIME_SLICE] < imarisTimeSlice]) - 1
    for i in range(slice_length):
        tmp.append(table[slice_length - i])
    return tmp


def rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells, table):
    if len(prevIterationSameLabelCells) != len(thisIterationSameLabelCells) and \
            [cell[CELL_ID], cell[CELL_LINE_NUMBER], cell[CELL_GENERATIONS], cell[CELL_TIME_SLICE]] not in \
            [
                [element[CELL_ID], element[CELL_LINE_NUMBER], element[CELL_GENERATIONS], element[CELL_TIME_SLICE]]
                for element in allPreviouslyLabeledCells(cell[CELL_TIME_SLICE], table)
            ]:
        return
    for labeledCell in prevIterationSameLabelCells:
        if not [labeledCell[CELL_LINE_NUMBER], labeledCell[CELL_LABEL]] in [[item[CELL_LINE_NUMBER], item[CELL_LABEL]] for item in thisIterationSameLabelCells]:
            cell[CELL_LABEL] = labeledCell[CELL_LABEL]
            cell[CELLS_LABELED_TOTAL] = labeledCell[CELLS_LABELED_TOTAL]
            break


XTTreeLabeled(0)
