#    <CustomTools>
#      <Menu>
#       <Item name="XTTreeLabeled" icon="Python">
#         <Command>Python3XT::XTTreeLabeled(%i)</Command>
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
CELL_DIVIDED = 5
CELL_TIME_SINCE_PREVIOUS_DIVISION = 6


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
    data = np.column_stack((statistics.mIds, statistics.mNames, statistics.mValues))
    generations = data[np.where(data[:, 1] == 'Generation')]
    timeIndex = data[np.where(data[:, 1] == 'Time Index')]
    timeSincePreviousDivision = data[np.where(data[:, 1] == 'Time Since Previous Division')]

    # Create table with CellId, TrackID, Generations, Time point
    table = prepareTableWithData(dataExtracted, generations, timeIndex, timeSincePreviousDivision)

    createUniqueLabels(table)

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
   :return [...[CELL_ID, CELL_LINE_NUMBER, CELL_GENERATIONS, CELL_TIME_SLICE, CELL_LABEL, CELL_DIVIDED]]
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
                table[j][CELL_DIVIDED] = 1
                break
            else:
                table[j][CELL_DIVIDED] = 0
    for j in range(len(table)):
        for i in range(len(timeSincePreviousDivision)):
            if int(timeSincePreviousDivision[i][CELL_ID]) == table[j][CELL_ID]:
                table[j][CELL_TIME_SINCE_PREVIOUS_DIVISION] = round(float(timeSincePreviousDivision[i][2]))
                break
            else:
                table[j][CELL_TIME_SINCE_PREVIOUS_DIVISION] = 0
    return table


def rewriteCells(currentSliceCells, previousIteration):
    assigned = []
    for currentSliceCell in currentSliceCells:
        for cell in previousIteration:
            if cell[CELL_LINE_NUMBER] == currentSliceCell[CELL_LINE_NUMBER] \
                    and cell[CELL_GENERATIONS] == currentSliceCell[CELL_GENERATIONS] \
                    and cell[CELL_LABEL] not in assigned \
                    and (
                    cell[CELL_TIME_SINCE_PREVIOUS_DIVISION] == (currentSliceCell[CELL_TIME_SINCE_PREVIOUS_DIVISION] - 1)
                    or (cell[CELL_GENERATIONS] == 0 and
                        cell[CELL_TIME_SINCE_PREVIOUS_DIVISION] == 0 and cell[CELL_TIME_SINCE_PREVIOUS_DIVISION] == currentSliceCell[
                            CELL_TIME_SINCE_PREVIOUS_DIVISION])
            ):
                currentSliceCell[CELL_LABEL] = cell[CELL_LABEL]
                assigned.append(cell[CELL_LABEL])
                break


def createUniqueLabels(table):
    """
   labels unique cells and adds number for cell during division
   :return void
    """
    previousIteration = []
    for timeSlice in range(int(max(table, key=lambda x: x[CELL_TIME_SLICE])[CELL_TIME_SLICE])):
        imarisTimeSlice = timeSlice + 1
        currentSliceCells = [cell for cell in table if int(cell[CELL_TIME_SLICE]) == imarisTimeSlice]
        currentSliceCells = sorted(currentSliceCells, key=lambda x: x[CELL_GENERATIONS])
        currentSliceCells = sorted(currentSliceCells, key=lambda x: x[CELL_TIME_SINCE_PREVIOUS_DIVISION])
        currentSliceCells = sorted(currentSliceCells, key=lambda x: x[CELL_DIVIDED])
        print(imarisTimeSlice)
        if imarisTimeSlice == 1:
            initializeLabels(currentSliceCells)
            previousIteration = currentSliceCells
            continue
        rewriteCells(currentSliceCells, previousIteration)
        for cell in currentSliceCells:
            if cell[CELL_LABEL] is not None:
                continue
            thisIterationSameLabelCells = [cellTmp for cellTmp in currentSliceCells if int(cellTmp[CELL_LINE_NUMBER]) == cell[CELL_LINE_NUMBER]]
            prevIterationSameLabelCells = [cellTmp for cellTmp in previousIteration if int(cellTmp[CELL_LINE_NUMBER]) == cell[CELL_LINE_NUMBER]]
            if cell[CELL_DIVIDED] == 0 and cell[CELL_LABEL] is None:
                rewriteFromAllLabeledCells(cell, table, thisIterationSameLabelCells, prevIterationSameLabelCells)
            if cell[CELL_LABEL] is not None:
                continue
            createNewLabel(cell, prevIterationSameLabelCells, thisIterationSameLabelCells)
        notLabeled = [e for e in currentSliceCells if e[CELL_LABEL] is None]
        if len(notLabeled) != 0:
            print(notLabeled)
            raise MissingLabelsException
        previousIteration = currentSliceCells


def createNewLabel(cell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    newLabel = getNewLabel(cell, prevIterationSameLabelCells, thisIterationSameLabelCells)
    cell[CELL_LABEL] = newLabel


def getNewLabel(cell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    usedLabels = [e[CELL_LABEL] for e in thisIterationSameLabelCells if e[CELL_LABEL]]
    usedLabels.sort()
    usedParents = list(set([e[:-2] for e in usedLabels]))
    usedParents.sort()
    parents = [p[CELL_LABEL] for p in prevIterationSameLabelCells if p[CELL_GENERATIONS] == (cell[CELL_GENERATIONS] - 1)]
    for parent in parents:
        if parent in usedLabels:
            continue
        candidate = parent + '_1'
        if candidate not in usedLabels:
            return candidate
        candidate = parent + '_2'
        if candidate not in usedLabels:
            return candidate

    print(cell)
    raise UnableToFindLabelException


def rewriteFromAllLabeledCells(cell, table, thisIterationSameLabelCells, previousIteration):
    rewriteCells(thisIterationSameLabelCells, previousIteration)
    notLabeled = [c for c in thisIterationSameLabelCells if c[CELL_LABEL] is None and c[CELL_DIVIDED] == 0]
    if not notLabeled:
        return
    allLabeledCells = list(set([c[CELL_LABEL] for c in table if
                                c[CELL_LABEL] and c[CELL_LINE_NUMBER] == cell[CELL_LINE_NUMBER] and c[CELL_GENERATIONS] == cell[CELL_GENERATIONS]]))
    usedLabels = [item[CELL_LABEL] for item in thisIterationSameLabelCells]
    for labeledCell in allLabeledCells:
        if labeledCell not in usedLabels:
            cell[CELL_LABEL] = labeledCell
            break


def initializeLabels(currentSliceCells):
    label = 1
    for cell in currentSliceCells:
        cell[CELL_LABEL] = str(label)
        label += 1


class UnableToFindLabelException(Exception):
    pass


class MissingLabelsException(Exception):
    pass


XTTreeLabeled(0)
