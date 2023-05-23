"""
This library contains useful function to process files and data commonly used in
OpenSim. It consists of Python implementation of functions that are already available
for MATLAB users, so when the authors were known their names have been added.
"""

import numpy as np
import scipy.signal as ss
from pathlib import Path

def readTRC(trcFile):
    """
    Read marker locations from .trc file. Optionally provides time, labels
    and units for marker locations.

    INPUTS
    :param trcFile: path to the .trc file to be read

    OUTPUTS
    :return markers: marker positions from the file
    :return times: vector of timestamps corresponding to marker trajectories
    :return labels: list of marker names
    :return units: string providing the units of measurements in which the data were
            saved

    USAGE: markers, times, labels, units = readTRC(trcFile)
    authors: Ajay Seth and Italo Belli
    """

    # extracting the main vectors from file
    markers_data = np.loadtxt(trcFile, delimiter='\t', skiprows=5)
    frameNums = markers_data[:, 0]
    times = markers_data[:, 1]
    markers = markers_data[:, 2:]

    # extracting labels and units of the .trc file
    with open(trcFile, 'r') as f:
        content = f.readlines()[1:4]
        f.close()
    fileProps = content[0].split()
    propValues = content[1].split()

    ixnm = fileProps.index('NumMarkers')
    ixu = fileProps.index('Units')
    numMarkers = propValues[ixnm]
    units = propValues[ixu]

    labels = content[2].split()
    labels = labels[2:]

    return markers, times, labels, units

def loadFilterCropArray(filepath, lowpassFreq, timeRange):
    """"
    This function allows to load data from a file (.sto or .mot) and to filter it
    with a 4th order Butterworth filter. The resulting data are cropped to return just
    those that are inside the required time range.

    INPUTS
    :param filepath: path to a STO or MOT file containing kinematic or force data
    :param lowpassFreq: lowpass frequency (Hz) used to process data with a 4th-order
           Butterworth filter
    :param timeRange: 1x2 array containing start and end time points where the data
           will be cropped.

    OUTPUTS
    :return array: filtered and cropped data array
    :return names: Data labels corresponding to the columns in 'array'
    :return time: Cropped time vector corresponding to the rows in 'array'

    authors: OpenSim Team Standford, Italo Belli (TU Delft)
    """

    # load data and convert it into Numpy arrays
    dataArray = np.loadtxt(filepath, delimiter='\t', skiprows=11)
    with open(filepath, 'r') as f:
        content = f.readlines()[0:11]
        f.close()
    names = content[10].split()
    ixtime = names.index('time')
    time = dataArray[:,ixtime]
    dataArray = np.delete(dataArray, ixtime, 1)
    names = np.delete(names, ixtime, 0)


    # Construct a 4th-order lowpass Butterworth filter
    # based on the 'lowpassFreq' argument
    timeStep = time[1]-time[0]
    sampleRate = 1/timeStep
    halfSampleRate = sampleRate/2
    cutoffFreq = lowpassFreq/halfSampleRate
    B,A = ss.butter(4, cutoffFreq, btype='low', output='ba')


    for i in range(0, dataArray.shape[1]-1):
        dataArray[:,i] = ss.filtfilt(B, A, dataArray[:,i], padtype='odd')

    # Extract subarrays based on 'timeRange' argument
    istart = (np.abs(time-timeRange[0])).argmin()
    iend = (np.abs(time-timeRange[1])).argmin()

    array = dataArray[istart:iend+1, :]
    time = time[istart:iend+1]

    return array, names, time


def writeMotStoData(data, optionType, columnNames, fileName):
    """
    This function allows to
    :param data: data to be written in the file
    :param optionType: set to 1 for a .mot output file, set to 2 for a .sto output file
    :param columnNames:  column headers in the resulting file
    :param fileName: name of the file that is generated
    :return:

    authors Adila Quesadilla Papaya and Italo Belli
    """
    nRows = data.shape[0]
    nColumns = data.shape[1]

    if optionType == 1:     # creating .mot file
        with open(fileName + '.mot', 'w+') as f:
            f.write(fileName + '\n')
            f.write('version=1 \n')
            f.write('nRows=' + str(nRows) + '\n')
            f.write('nColumns=' + str(nColumns) + '\n')
            f.write('inDegrees=yes \n')
            f.write('endheader\n')
            f.write('time\t')

            for i in range(0, nColumns-1):
                f.write(columnNames[i] + ' \t ')

            f.write('\n')

            for i in range(0, nRows):
                for j in range(0, nColumns):
                    f.write('\t ' + str("{:.6f}".format(data[i, j])))
                f.write('\n')

            f.write('\n')
            f.close()

    if optionType == 2:     # creating .sto file
        with open(fileName + '.sto', 'w+') as f:
            f.write(fileName + '\n')
            f.write('version=1 \n')
            f.write('nRows=' + str(nRows) + '\n')
            f.write('nColumns=' + str(nColumns) + '\n')
            f.write('inDegrees=yes \n')
            f.write('endheader\n')
            f.write('time\t')

            for i in range(0, nColumns-1):
                f.write(columnNames[i] + ' \t ')

            f.write('\n')

            for i in range(0, nRows):
                for j in range(0, nColumns):
                    f.write('\t ' + str("{:.6f}".format(data[i, j])))
                f.write('\n')

            f.write('\n')
            f.close()
