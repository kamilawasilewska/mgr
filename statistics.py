import os

import matplotlib.pyplot as plot
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

# path = "C:\\Users\\Kamila\\Documents\\005 FITC_C_A-80_Statistics"
path = "C:\\Users\\Kamila\\Documents\\005 FITC_C_A_all_slice_Statistics"
arr: list = os.listdir(path)


def get_file_names_with_strings(str_list: list):
    full_list = arr
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list


filename = get_file_names_with_strings(['Intensity_Center_Ch=2'])
treeLabeled: list = get_file_names_with_strings(['Tree_Labeled'])

filename = str(path) + '\\' + str(filename[0])
data: pd.DataFrame = pd.read_table(filename, header=2, sep=",")

treeLabeledData: str = str(path) + '\\' + str(treeLabeled[0])
treeLabeledData: pd.DataFrame = pd.read_table(treeLabeledData, header=2, sep=",")

imarisParsedData: np.ndarray = np.column_stack(
    (treeLabeledData['ID'], data['Intensity Center'], data['Time'], treeLabeledData['TreeLabeled']))

cellLabels = np.unique(list(imarisParsedData[:, 3]))
cellsWithIntensityByTime: dict = {}
cellGenerations: dict = {}
allCells: dict = {}
for cellLabel in cellLabels:
    cellsWithIntensityByTime[cellLabel] = {entry[0]: [entry[2], entry[1]] for entry in imarisParsedData if
                                           entry[3] == cellLabel}
    allCells[cellLabel] = [[cell[2], cell[1]] for cell in imarisParsedData if cell[3] == cellLabel]
    labelAsArray: imarisParsedData = cellLabel.split('_')
    if len(labelAsArray) == 1 or labelAsArray[-1] != '1':
        continue
    itemIndex: np.ndarray = np.where(imarisParsedData == cellLabel)[0][0]
    if len(cellGenerations) == 0:
        parentCandidates = [e for e in imarisParsedData if e[2] == imarisParsedData[itemIndex][2]]
        parent: list = [parent for parent in parentCandidates if
                        parent[-1] == imarisParsedData[itemIndex][-1].split('_')[0]]
        idx = itemIndex
        while not parent:
            idx -= 1
            parentCandidates = [e for e in imarisParsedData if e[2] == imarisParsedData[idx][2]]
            parent = [parent for parent in parentCandidates if
                      parent[-1] == imarisParsedData[itemIndex][-1].split('_')[0]]
            if parent:
                break

        cellGenerations[parent[0][-1]] = {parent[0][0]: [parent[0][2], parent[0][1]]}
    tmp = imarisParsedData[itemIndex]
    cellGenerations[cellLabel] = {tmp[0]: [tmp[2], tmp[1]]}

dividingCellLabels: list = list(cellGenerations.keys())
allValues: list = []
for dividingCellLabel in dividingCellLabels:
    tmp = [[dataRow[2], dataRow[1]] for dataRow in imarisParsedData if dataRow[3] == dividingCellLabel]
    timeSlice = [item[0] for item in allValues]
    for i in tmp:
        if not i[0] in timeSlice:
            allValues.append(i)

x: list = []
y: list = []
cellID: tuple
unused: tuple
cellID, unused = zip(*list(cellGenerations.items()))
for element in [[(k, v) for k, v in value.items()] for value in cellGenerations.values()]:
    x.append(element[0][1][0])
    y.append(element[0][1][1])

fix: Figure
fix, ax = plot.subplots()
ax.plot(x, y)
ax.set_xlabel('czas (slice)')
ax.set_ylabel('wartość świecenia')
ax.set_title('Wykres zmian świecenia względem pokolenia')
for i in range(len(cellID)):
    plot.text(x[i], y[i], cellID[i])
plot.show()

fix, ax = plot.subplots()
ax.plot([item[0] for item in allValues], [item[1] for item in allValues])
ax.set_xlabel('czas (slice)')
ax.set_ylabel('wartość świecenia')
ax.set_title('Wykres zmian świecenia w czasie')
for i in range(len(cellID)):
    plot.text(x[i], y[i], cellID[i])
plot.show()

legend: list = []
x: list = []
y: list = []
ids: list = []
cellIds, unused = zip(*list(allCells.items()))
idx = 0
for data in allCells.items():
    plot.plot([item[0] for item in data[1]], [item[1] for item in data[1]])
    #legend.append('Cell: ' + data[0])
    if not data[0] in ids:
        ids.append(cellIds[idx])
        x.append(unused[idx][0][0])
        y.append(unused[idx][0][1])
        idx += 1
for i in range(len(ids)):
    plot.scatter(x[i], y[i])
    #plot.text(x[i], y[i], ids[i])

plot.title('Wykres zmian świecenia w czasie (wszystkie komórki)')
plot.ylabel('wartość świecenia')
plot.xlabel('czas (slice)')
plot.legend(legend)
plot.show()
