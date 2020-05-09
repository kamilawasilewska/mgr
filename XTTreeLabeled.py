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
        if imarisTimeSlice == 1:
            label = 1
            for cell in currentSliceCells:
                cell[4] = str(label)
                cell[5] = 1
                label += 1
        else:
            for cell in currentSliceCells:
                thisIterationSameLabelCells = [l for l in currentSliceCells if int(l[1]) == cell[1]]
                prevIterationSameLabelCells = [l for l in previousIteration if int(l[1]) == cell[1]]
                if len(prevIterationSameLabelCells) < len(thisIterationSameLabelCells):
                    rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells)
                    if cell[4] is None:
                        try:
                            label
                        except NameError:
                            label = 1
                        try:
                            cell[4] = prevIterationSameLabelCells[0][4] + '_' + str(label)
                        except TypeError:
                            tmp = []
                            slice_length = len([l for l in table if l[3] < imarisTimeSlice]) - 1
                            for i in range(slice_length):
                                tmp.append(table[slice_length - i])
                            for labeledCell in tmp:
                                if not labeledCell[5] in [item[-1] for item in thisIterationSameLabelCells]:
                                    cell[4] = prevIterationSameLabelCells[0][4] + '_' + label
                                    break
                        del tmp, slice_length
                        cell[5] = label
                        label += 1
                elif len(prevIterationSameLabelCells) >= len(thisIterationSameLabelCells):
                    rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells)
            try:
                del label
            except NameError:
                print('no label to delete')
        previousIteration = currentSliceCells


def rewriteLabeledCells(cell, prevIterationSameLabelCells, thisIterationSameLabelCells):
    for labeledCell in prevIterationSameLabelCells:
        if not [labeledCell[1], labeledCell[5]] in [[item[1], item[-1]] for item in thisIterationSameLabelCells]:
            cell[4] = labeledCell[4]
            cell[5] = labeledCell[5]
            break
