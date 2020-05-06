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
from collections import Counter
from tkinter import *


def XTTreeLabeled(imarisId):
    # Create an ImarisLib object
    imarisLib = ImarisLib.ImarisLib()

    # Get an imaris object with id aImarisId
    imarisId = 0
    imarisApp = imarisLib.GetApplication(imarisId)

    # Get the currently loaded dataset
    currentImage = imarisApp.GetImage(0)

    # Get Cell from Imaris
    # cells = imarisApp.GetFactory().ToCells(imarisApp.GetSurpassSelection())# Works only for Cell mode

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

    #     mFactorsCells = statistics.mFactors
    #     mFactorNameCells = statistics.mFactorNames
    #     mIdsCells = statistics.mIds #This value has tha same result like function object.GetIds()
    #     mNamesCells = statistics.mNames
    #     mValuesCells = statistics.mValues

    data = np.column_stack((statistics.mIds, statistics.mNames, statistics.mValues))
    # uniqueSpecificName = np.transpose(np.unique(np.array(statistics.mNames)))
    generations = data[np.where(data[:, 1] == 'Generation')]
    timeIndex = data[np.where(data[:, 1] == 'Time Index')]

    dataExtracted = np.array([[15, 100015],
                              [24, 100024],
                              [64, 100015],
                              [70, 100024],
                              [78, 100015],
                              [80, 100024],
                              [90, 100024],
                              [92, 100015],
                              [95, 100024],
                              [96, 100024]
                              ])

    generations = np.array([['15', 'ok', 0],
                            ['24', 'ok', 0],
                            ['64', 'ok', 0],
                            ['70', 'ok', 0],
                            ['78', 'ok', 1],
                            ['80', 'ok', 1],
                            ['90', 'ok', 1],
                            ['92', 'ok', 1],
                            ['95', 'ok', 0],
                            ['96', 'ok', 1]
                            ])

    timeIndex = np.array([['15', 'ok', 1],
                          ['24', 'ok', 1],
                          ['64', 'ok', 2],
                          ['70', 'ok', 2],
                          ['78', 'ok', 3],
                          ['80', 'ok', 3],
                          ['90', 'ok', 3],
                          ['92', 'ok', 4],
                          ['95', 'ok', 4],
                          ['96', 'ok', 4]
                          ])
    # Create table with CellId, TrackID, Generations, Time point
    table = prepareTableWithData(dataExtracted, generations, timeIndex)
    res = np.array(table)

    createUniqueLabels(table)

    res = np.array(table)


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
        if dataExtracted[i, 0] == int(generations[i, 0]) and dataExtracted[i, 0] == int(timeIndex[i, 0]):
            table.append([dataExtracted[i, 0], dataExtracted[i, 1], int(float(generations[i, 2])), int(float(timeIndex[i, 2])), None, None])
    return table


def createUniqueLabels(table):
    """
   labels unique cells and adds number for cell during division
   :return void
    """
    for timeSlice in range(int(max(table, key=lambda x: x[3])[3])):
        imarisTimeSlice = timeSlice + 1
        currentSliceCells = [l for l in table if int(l[3]) == imarisTimeSlice]
        if imarisTimeSlice == 1:
            label = 1
            for cell in currentSliceCells:
                cell[4] = str(label)
                cell[5] = 1
                label += 1
        else:
            for cell in currentSliceCells:
                hasCellDivided = len([l for l in currentSliceCells if int(l[1]) == cell[1]]) > 1
                prevIterationCells = [l for l in previousIteration if int(l[1]) == cell[1]]
                cells = [l for l in currentSliceCells if int(l[1]) == cell[1]]
                if len(prevIterationCells) < len(cells):
                     #Przypadek gdy na poprzednim slajsie jest mniejsza ilość komórek niż na obecnym. (Podział komórek)
                    try:
                        label
                    except NameError:
                        label = 1
                    else:
                        label += 1
                    try:
                        cell[4] = prevIterationCells[0][4] + '_' + str(label)
                    except IndexError:
                        cell[4] = str(label)
                    except TypeError:
                        tmp = []
                        for i in range(timeSlice):
                            tmp.append(table[timeSlice - i])
                        for labeledCell in tmp:
                            if labeledCell[1] == cell[1] and not labeledCell[5] in [item[-1] for item in
                                                                                    currentSliceCells]:
                                cell[4] = labeledCell[4]
                                break
                    cell[5] = label
                elif len(prevIterationCells) > len(cells) and not hasCellDivided:
                     #Ilość obecnych komórek jest mniejsza niż na poprzednim slajsie. I nie nastąpił ich podział.
                     #wycinek glównej table timeSlice do 0
                     #i szukać wiersza w którym table[row,1] == cell[1]
                    #w 5 kolumnie musi być unikalna wartośc dla tego wycinku czasowego
                    tmp = []
                    for i in range(timeSlice):
                        tmp.append(table[timeSlice - i])
                    for labeledCell in tmp:
                        if labeledCell[1] == cell[1] and not labeledCell[5] in [item[-1] for item in currentSliceCells]:
                            cell[4] = labeledCell[4]
                            cell[5] = labeledCell[5]
                            break
                elif len(prevIterationCells) == len(cells) and hasCellDivided:
                    #Nastąpił podział komórek, ale któraś musiała umrzeć.
                    try:
                        label
                    except NameError:
                        label = 1
                    else:
                        label += 1
                    prevCell = [l for l in prevIterationCells if int(l[5]) == label][0]
                    cell[4] = prevCell[4]
                    cell[5] = prevCell[5]
                elif hasCellDivided:
                    #Tylko nastąpił podział komórek, nie sprawdzamy liczebności poprzedniej i tej.
                    try:
                        label
                    except NameError:
                        label = 1
                    else:
                        label += 1
                    try:
                        cell[4] = prevIterationCells[0][4] + '_' + str(label)
                    except TypeError:
                        tmp = []
                        for i in range(timeSlice):
                            tmp.append(table[timeSlice - i])
                        for labeledCell in tmp:
                            if labeledCell[1] == cell[1] and not labeledCell[5] in [item[-1] for item in
                                                                                    currentSliceCells]:
                                cell[4] = labeledCell[4]
                                break
                    cell[5] = label
                else:
                    #Stan komórek bez zmian
                    cell[4] = prevIterationCells[0][4]
                    cell[5] = prevIterationCells[0][5]
            try:
                del label
            except NameError:
                pass
            except IndexError:
                pass
        previousIteration = currentSliceCells