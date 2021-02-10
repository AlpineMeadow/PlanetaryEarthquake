#!/usr/bin/env python3


# Gather our code in a main() function
def main():
    import os.path
    import PlanetaryEarthquakeFunctions as pf
    import numpy as np
    import math
    import matplotlib.pyplot as plt
#    matplotlib inline
#    from mpl_toolkits.basemap import Basemap, cm


    #Set up the path and file name for the earthquake data.
    filepath = '/home/jdw/UM2019Autumn/M561/Project/'
    fname = 'SolarSystemAndEarthquakes.csv'
    earthquakeFilename = filepath + fname

    #Get the dictionary that will hold the earthquake data.
    eqdata = pef.getEarthquakeData(earthquakeFilename)

    #Get the ephemeris data and the force data.
    #Check to see if the file holding the data has already
    #been produced.
    savefilename = 'ProjectSaveData.npy'
    savefile = filepath + savefilename
    if os.path.isfile(savefile) :
        data = np.load(savefile)    
    else :
        data = pef.getData(eqdata)
        np.save(savefile, data)


    #Make histograms of the earthquake data.
    pef.eqHistograms(data)
    
    #Make plots of force versus latitude for each planet.
    pef.polarplot(data)

    #Make plots of earthquake location and planet location
    #for various times.
    pef.EqPlanetDistance(data)
    
    #Make plots of force versus time for each earthquake.
    pef.forcePlot(data)

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
  main()

