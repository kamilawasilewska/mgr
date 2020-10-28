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

    division = []
    for i in range(len(table)):
        if table[i][CELL_TIME_SINCE_PREVIOUS_DIVISION] == 2:
            division.append(table[i])

    createLabels(table)
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


def createLabelForPreviouslyDividedCell(currentSliceCells, table):
    previouslyDividedNotLabeled = [cell for cell in currentSliceCells if cell[CELL_LABEL] is None and cell[CELL_DIVIDED] == 0]
    for missingLabelCell in previouslyDividedNotLabeled:
        tmp = [cell for cell in table if
               cell[CELL_LINE_NUMBER] == missingLabelCell[CELL_LINE_NUMBER] and cell[CELL_GENERATIONS] == missingLabelCell[CELL_GENERATIONS] - 1]
        tmp = [cell[CELL_LABEL] for cell in tmp if cell[CELL_LABEL] is not None]
        parentCandidates = list(set(tmp))
        parentCandidates.sort()
        labelsTaken = [
            cell[CELL_LABEL] for cell in currentSliceCells if
            cell[CELL_LABEL] is not None and cell[CELL_LINE_NUMBER] == missingLabelCell[CELL_LINE_NUMBER] and cell[CELL_GENERATIONS] == missingLabelCell[
                CELL_GENERATIONS]
        ]

        tmp = []
        for parentLabel in labelsTaken:
            tmp.append(parentLabel[:-2])

        realCandidates = []
        for item in tmp:
            if tmp.count(item) == 1:
                realCandidates.append(item)
                realCandidates.append(item)

        label = None
        for newLabel in realCandidates:
            if label is None or label > 2:
                label = 1
            newLabel += '_'
            newLabel += str(label)
            if newLabel in labelsTaken:
                label += 1
                continue
            missingLabelCell[CELL_LABEL] = newLabel
            break


# 2122

def createLabels(table):
    usedLabels = []
    for timeSlice in range(int(max(table, key=lambda timeWindow: timeWindow[CELL_TIME_SLICE])[CELL_TIME_SLICE])):
        imarisTimeSlice = timeSlice + 1
        print(imarisTimeSlice)
        label = None
        currentSliceCells = [cell for cell in table if int(cell[CELL_TIME_SLICE]) == imarisTimeSlice]
        currentSliceCells = sorted(currentSliceCells, key=lambda x: x[CELL_GENERATIONS])
        if imarisTimeSlice == 1:
            initializeLabels(currentSliceCells)
            continue

        rewriteLabelBasingOnTimeSincePreviousDivision(currentSliceCells, table)
        rewriteLabelForCells(currentSliceCells, table)
        createLabelForPreviouslyDividedCell(currentSliceCells, table)
        createNewLabels(currentSliceCells, label, table, usedLabels)


def createNewLabels(currentSliceCells, label, table, usedLabels):
    newCells = [newCell for newCell in currentSliceCells if newCell[CELL_DIVIDED] == 1]
    tmp = [cell[CELL_LABEL] for cell in currentSliceCells if cell[CELL_LABEL] is not None]
    usedLabels = list(set(tmp + usedLabels))
    usedLabels.sort()
    for newCell in newCells:
        parentLabel = findParentLabel(newCell, usedLabels, currentSliceCells, table)
        label = createLabel(newCell, parentLabel, label)
        usedLabels.append(parentLabel)


def rewriteLabelBasingOnTimeSincePreviousDivision(currentSliceCells, table):
    # todo tutaj w 113 powinien poprawnie przepisac komórki a tego nie robi
    # powinien przepisać wszystkie poza tymi ktore sie dziela
    notDividedCells = [cell for cell in currentSliceCells if cell[CELL_LABEL] is None and cell[CELL_DIVIDED] == 0]
    notDividedCells.sort(key=lambda cell: cell[CELL_TIME_SINCE_PREVIOUS_DIVISION])
    previousIterationCells = [prevCell for prevCell in table if prevCell[CELL_TIME_SLICE] == notDividedCells[0][CELL_TIME_SLICE] - 1]
    previousIterationCells.sort(key=lambda cell: cell[CELL_TIME_SINCE_PREVIOUS_DIVISION])
    usedLabels = []
    for notDividedCell in notDividedCells:
        for prevCell in previousIterationCells:
            if (
                    # cell is not dividing just rewrite
                    notDividedCell[CELL_LINE_NUMBER] == prevCell[CELL_LINE_NUMBER]
                    and notDividedCell[CELL_GENERATIONS] == prevCell[CELL_GENERATIONS]
                    and notDividedCell[CELL_TIME_SINCE_PREVIOUS_DIVISION] == prevCell[CELL_TIME_SINCE_PREVIOUS_DIVISION]
                    and notDividedCell[CELL_TIME_SINCE_PREVIOUS_DIVISION] == 0
                    and prevCell[CELL_LABEL] not in usedLabels
            ) \
                    or (
                    # cell present in previous iteration also rewrite
                    notDividedCell[CELL_LINE_NUMBER] == prevCell[CELL_LINE_NUMBER]
                    and notDividedCell[CELL_GENERATIONS] == prevCell[CELL_GENERATIONS]
                    and (notDividedCell[CELL_TIME_SINCE_PREVIOUS_DIVISION] - 1) == prevCell[CELL_TIME_SINCE_PREVIOUS_DIVISION]
                    and prevCell[CELL_LABEL] not in usedLabels
            ):
                notDividedCell[CELL_LABEL] = prevCell[CELL_LABEL]
                usedLabels.append(prevCell[CELL_LABEL])
                break


def initializeLabels(currentSliceCells):
    label = 1
    for cell in currentSliceCells:
        cell[CELL_LABEL] = str(label)
        label += 1


def rewriteLabelForCells(currentSliceCells, table):
    notDividedCellsWithoutLabel = [tmp for tmp in currentSliceCells if tmp[CELL_DIVIDED] == 0 and tmp[CELL_LABEL] is None]
    if not notDividedCellsWithoutLabel:
        return
    labeledCells = [labeledCell for labeledCell in table if int(labeledCell[CELL_TIME_SLICE]) < notDividedCellsWithoutLabel[0][CELL_TIME_SLICE]]
    for cellToRewrite in notDividedCellsWithoutLabel:
        for labeledCell in labeledCells:
            if cellToRewrite[CELL_GENERATIONS] == labeledCell[CELL_GENERATIONS] \
                    and cellToRewrite[CELL_LINE_NUMBER] == labeledCell[CELL_LINE_NUMBER] \
                    and not labeledCell[CELL_LABEL] in [item[CELL_LABEL] for item in notDividedCellsWithoutLabel] \
                    and not labeledCell[CELL_LABEL] in [currentSliceCell[CELL_LABEL] for currentSliceCell in currentSliceCells if
                                                        currentSliceCell[CELL_LABEL] is not None]:
                cellToRewrite[CELL_LABEL] = labeledCell[CELL_LABEL]
                break


def findParentLabel(cell, usedLabels, rewrittenCells, table):
    # tu jest ta kurwa zjebana na 113 iteracji nadaje 2x ten sam label
    alreadyUsedLabels = [tmp[CELL_LABEL] for tmp in rewrittenCells if tmp[CELL_LABEL] is not None and cell[CELL_LINE_NUMBER] == tmp[CELL_LINE_NUMBER]]
    alreadyUsedLabels.sort()
    labels = [
        labeledCell for labeledCell in table if
        int(labeledCell[CELL_TIME_SLICE]) < cell[CELL_TIME_SLICE] \
        and labeledCell[CELL_LINE_NUMBER] == cell[CELL_LINE_NUMBER] \
        and labeledCell[CELL_GENERATIONS] == (cell[CELL_GENERATIONS] - 1)
        and labeledCell[
            CELL_LABEL] not in alreadyUsedLabels
    ]
    labels = list(set([cell[CELL_LABEL] for cell in labels]))
    labels.sort()
    tmp = []
    for candidate in labels:
        if candidate is not None and len(candidate.split('_')) == cell[CELL_GENERATIONS]:
            tmp.append(candidate)
    tmp.sort()
    for labeledCell in reversed(tmp):
        if usedLabels.count(labeledCell) == 2:
            continue
        else:
            return labeledCell

    a = 1


def createLabel(cell, parentLabel, label):
    if label is None or label > 2:
        label = 1
    parentLabel += '_'
    parentLabel += str(label)
    cell[CELL_LABEL] = parentLabel
    label += 1
    return label


XTTreeLabeled(0)
