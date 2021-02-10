#A module containing functions needed for the M561 project.



#A function that will get the earthquake data.
def getEarthquakeData(filename) :
    import spiceypy as spice

    #Load the meta kernel file for the conversion from UTC to ET.
    metakernel = '/home/jdw/UM2019Autumn/M561/Project/kernels/M561ProjectKernels.txt'
    spice.furnsh( metakernel )

    #Create a dictionary to hold the earthquake data.
    eqdict = {}

    #Open the Earthquake file and get the time, location and magnitude.
    with open(filename , 'r') as f :
        header_line = f.readline()
        
        #The data I want is in the following positions:
        #data[0] holds the time in UTC.
        #data[1] holds the earthquake latitude (units of degrees).
        #data[2] holds the earthquake longitude (units of degrees).
        #data[3] holds the earthquake magnitude.
        for string in f :

            data = string.split(",")
            
            #Convert the time information into a string to be input into
            #the spiceypy function.
            t_str = str(data[0][0 : 4] + ' ' + data[0][5 : 7] + ' ' + data[0][8 : 10] + ' ' +
                     data[0][11 : 13] + ':' + data[0][14 : 16] + ':' + data[0][17 : 23])
            
            #Call the spiceypy function that converts a normal time string
            #into an ephemeris time string.
            et = spice.str2et ( t_str )

            #Fill the dictionary.
            eqdict[et] = (data[3], data[1], data[2])

    #Close the file.
    f.close

    #Return the dictionary.
    return eqdict

############################################################

#A function that will generate the distance from Earth to planet
#for the given earthquake times.
def getData(eqdata) :
    import spiceypy as spice
    import numpy as np
    import constants as con
    import math
    
    #Load the meta kernel file for the determination of the position of the planets.
    metakernel = '/home/jdw/UM2019Autumn/M561/Project/kernels/M561ProjectKernels.txt'
    spice.furnsh( metakernel )

    #Find the number of earthquakes to be analyzed.
    nquake = len(eqdata)

    #Create a numpy array with nquake rows and 49 columns.
    data = np.zeros( (nquake, 284), dtype = 'float64' )
    
    #Get the ephemeris times from the earthquake dictionary keys.  
    et = list(eqdata.keys())
    
    #Get the earthquake data from the earthquake dictionary values.
    eqdata = list(eqdata.values())

    #Get the positions
    #First create a vector of planet names to be fed to the spice
    #function.
    Planets = ['Sun','Mercury_Barycenter', 'Venus_Barycenter', 'Moon',
               'Mars_Barycenter', 'Jupiter_Barycenter', 'Saturn_Barycenter',
               'Uranus_Barycenter', 'Neptune_Barycenter', 'Pluto_Barycenter']
    
    #Create a time offset index.  I want to look not only at the exact earthquake
    #time but also times that are a given number of seconds before and after
    #the earthquake.
    k = [3000, 2400, 1800, 1200, 600, 0, -600]
    
    for i in range(nquake) :
        #First fill in the earthquake data.
        data[i][0] = et[i]
        data[i][1] = eqdata[i][0]  #Earthquake Magnitude
        data[i][2] = eqdata[i][1]  #Earthquake Latitude (degrees)
        data[i][3] = eqdata[i][2]  #Earthquake Longitude (degrees)

        #Now loop through the planets to get the distances.
        for p in range(len(Planets)) :
            pindex = p*28 + 4

            #Now loop through the seven times around the earthquake.
            for t in range(7) :
                
                #Create a time value.  The earthquakes are given in the et
                #vector and I want to look at times before those earthquake
                #times.
                time = et[i] - k[t]
                
                #Determine the position of a given planet with respect to Earth
                #for a given et in the J2000 coordinate system.
                pos, ltime = spice.spkpos( Planets[p], time, 'J2000', 'NONE', 'Earth', )

                #Calculate the distance in meters. Vnorm returns a distance in kilometers
                #so we have to convert.  1km = 1000.0 meters.
                dist = spice.vnorm(pos)*1000.0

                #Return the range, longitude and latitude of the position vector
                #pointing from Earth to the planet.  These are output in radians.
                ran, lon, lat = spice.reclat(pos)
                
                #Convert the radians to degrees.
                longitude = math.degrees(lon)
                latitude = math.degrees(lat)

                #Get the gravitational force for the times preceeding and
                #following the earthquake.
                F = getGravForce(p, dist)
                
                #Fill the data array.
                data[i][pindex + t*4] = dist
                data[i][pindex + t*4 + 1] = latitude
                data[i][pindex + t*4 + 2] = longitude
                data[i][pindex + t*4 + 3] = F
    
    return data


###########################################################


#A function that will determine the force on Earth by the planets
#for a given earthquake time.
#I will determine if there needs to be a tidal force by calculating the
#forces from the edge of the planet and comparing those forces to
#that force from just the center of the planet.
def getGravForce(p, dist) :
    import constants as con
    import spiceypy as spice
    
    Planets = ['Sun','Mercury', 'Venus', 'Moon', 'Mars', 'Jupiter', 'Saturn',
               'Uranus', 'Neptune', 'Pluto']

    GMe = con.G*con.PlanetaryData['Earth'][0]
    d = dist

    #Calculate the force.  Assume that there is no tidal force.
    Force = GMe*con.PlanetaryData[Planets[p]][0]/(d*d)
    
    return Force

###########################################################


def eqHistograms(data) :
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import spiceypy as spice
    from matplotlib.backends.backend_pdf import PdfPages

    #Turn the ephemeris time into Universal Time.
    calet = spice.etcal( data[0][0] )

    #Get the data for the Earthquake magnitude, Earthquake latitude
    #and Earthquake longitude.
    mag = data[ : , 1].copy()
    lat = data[ : , 2].copy()
    lon = data[ : , 3].copy()
    
    #Plot a histogram of the magnitude of the earthquakes.
    outfilepath = '/home/jdw/UM2019Autumn/M561/Project/'
    outfile = outfilepath + 'Plots/EQMagHist.pdf'

    n, bins, patches = plt.hist(mag, 25, density = 1, facecolor='g')
    plt.xlabel('Magnitude')
    plt.ylabel('Probability')
    plt.title('Histogram of Earthquake Magnitude')
    plt.grid(True)
        
    #Save the plot to a file.
    plt.savefig(outfile)

    #Plot a histogram of the latitude of the earthquakes.
    outfile = outfilepath + 'Plots/EQLatHist.pdf'
    n, bins, patches = plt.hist(lat, 50, density=1, facecolor='g', alpha=0.75)
    plt.xlabel('Latitude')
    plt.ylabel('Probability')
    plt.title('Histogram of Earthquake Latitude')
    plt.axis([-90, 90, 0, 0.03])
    plt.grid(True)
    
    #Save the plot to a file.
    plt.savefig(outfile)

    
    #Plot a histogram of the longitude of the earthquakes.
    outfile = outfilepath + 'Plots/EQLongHist.pdf'
    n, bins, patches = plt.hist(lon, 50, density=1, facecolor='g', alpha=0.75)
    plt.xlabel('Longitude')
    plt.ylabel('Probability')
    plt.title('Histogram of Earthquake Longitude')
    plt.axis([-180, 180, 0, 0.02])
    plt.grid(True)

    #Save the plot to a file.
    plt.savefig(outfile)

##############################################################

#Make polar plots of force versus latitude for each planet.
def polarplot(data) :
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import cartopy.crs as ccrs
    import spiceypy as spice
    from matplotlib.backends.backend_pdf import PdfPages

    #Get the data for the Earthquake magnitude, Earthquake latitude
    #and Earthquake longitude.
    mag = data[ : , 1].copy()
    lat = data[ : , 2].copy()
    lon = data[ : , 3].copy()
    
    #Plot a histogram of the magnitude of the earthquakes.
    outfilepath = '/home/jdw/UM2019Autumn/M561/Project/'
    outfile = outfilepath + 'Plots/GlobalEqPlanetLocation.png'

    #Find the number of rows in the data matrix.
    nRows = data.shape[0]
    print(nRows)
    #Allocate an array of distances for each planet for each earthquake.
    LatPlanet = np.zeros((nRows, 10), dtype = 'float64')
    LonPlanet = np.zeros((nRows, 10), dtype = 'float64')

    #Allocate vectors of earthquake latitudes and longitudes.
    EQLat = np.zeros((nRows), dtype = 'float64')
    EQLon = np.zeros((nRows), dtype = 'float64')
    calet = [0]*nRows
    
    for i in range(nRows) :
        #Turn the ephemeris time into Universal Time.
        calet[i] = spice.etcal( data[i][0] )
        EQLat[i] = data[i, 2] #given in degrees.
        EQLon[i] = data[i, 3] #given in degrees.
        for j in range(10) :
            index = 25 + j*28
            LatPlanet[i, j] = data[i, index] #given in degrees.
            LonPlanet[i, j] = data[i, index + 1] #given in degrees.

    #Choose some times to look at earthquakes.
    time_vec = [300, 900, 1500, 2100, 2700]
    
    plt.figure()
    #        plt.gcf()
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines()
    ax.stock_img()

    titlestr = (r"Earthquake and Planet Locations for " +
                "\n" + "M = " + str(data[300, 1]) +
                " Earthquake on " + "\n" + calet[300])

    #Plot the location of the earthquake.
    plt.scatter(EQLon[300], EQLat[300], color='red', marker='o', transform=ccrs.Geodetic())
    plt.title(titlestr)

    #Plot the locations of the planets.
    plt.scatter(LatPlanet[300, -10:], LonPlanet[300, -10:], color='blue', marker='o',
                transform=ccrs.Geodetic())

    #Enter some text
    plt.text(EQLon[300] - 1, EQLat[300] - 10, 'Earthquake', horizontalalignment='right',
             color = 'red', transform=ccrs.Geodetic())

    plt.text(LatPlanet[300, 0] - 6, LonPlanet[300, 0] - 12, 'Sun',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 1] + 2, LonPlanet[300, 1] + 4, 'Mercury',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 2] + 1, LonPlanet[300, 2] - 6, 'Venus',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 3] + 3, LonPlanet[300, 3] - 12, 'Moon',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 4] + 3, LonPlanet[300, 4] + 4, 'Mars',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 5] - 3, LonPlanet[300, 5] - 12, 'Jupiter',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 6] + 2, LonPlanet[300, 6] + 3, 'Saturn',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 7] + 3, LonPlanet[300, 7] - 12, 'Uranus',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 8] + 3, LonPlanet[300, 8] - 12, 'Neptune',
             horizontalalignment='left', transform=ccrs.Geodetic())
    
    plt.text(LatPlanet[300, 9] + 3, LonPlanet[300, 9] - 12, 'Pluto',
             horizontalalignment='left', transform=ccrs.Geodetic())

    #Save the plot to a file.
    plt.savefig(outfile)





#############################################################    

def function_plot(X, color, planet, pp):
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    import numpy as np
    
    # Fit a normal distribution to the data:
    mu, sigma = norm.fit(X)
    
    #Create a text string containing the values of the mean and standard
    #deviation.
    textstr = '\n'.join((r'$\mu=%.2f$ km' % (mu, ), r'$\sigma=%.2f$ km' % (sigma, )))

    #Set the number of bins
    numbins = 100

    #Generate the figure instance.
    plt.figure()
    ax = plt.axes()

    #Clear the figure.
    plt.clf()

    #Plot the histogram.
    plt.hist(X, color = color, label = planet, bins = numbins, density = True)

    #Find the minimum and maximum x and y values.
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    #Create a vector of evenly space points spanning from xmin to xmax.
    xspace = np.linspace(xmin, xmax, numbins)
    
    #Create the probability density function for the given mu and sigma.
    p = norm.pdf(xspace, mu, sigma)

    #Plot the normal curve.
    plt.plot(xspace, p, 'k', linewidth=1)
    plt.title('Histogram of Distance from Earthquake to ' + planet)
    plt.ylabel('Relative Frequency')
    plt.xlabel('Distance (km)')
    plt.yticks(fontsize=6)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

    # place a text box in upper left in axes coords
    plt.text(0,0.8*ymax, textstr, fontsize=14, bbox=props)

#    pp.savefig(plt.gcf())

    outfilepath = '/home/jdw/UM2019Autumn/M561/Project/'
    outfile = outfilepath + 'Plots/EqPlanetDistance.png'
    plt.savefig(outfile)

#############################################################
    
#Make plots of earthquake location and planet location.
def EqPlanetDistance(data) :
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import constants as con
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.backends.backend_pdf import PdfPages
    
    #Get the radius of Earth.
    Re = con.PlanetaryData['Earth'][1]

    #Create a string holding the location of project.
    outfilepath = '/home/jdw/UM2019Autumn/M561/Project/'

    #Create a file name for the output plot.
    outfile = outfilepath + 'Plots/EqPlanetDistance.pdf'
    
    #Find the number of rows in the data matrix.
    nRows = data.shape[0]

    #Allocate an array of distances for each planet for each earthquake.
    distanceEqPlanet = np.zeros((nRows, 10), dtype = 'float64')
    
    #Now loop though the data to get the distances..
    for i in range(nRows) :
        #Set up some planetary latitude and longitude vectors.
        EqLat = math.radians(data[i,2])
        EqLon = math.radians(data[i,3])
        
        #Loop through the planets to get the latitude and longitude
        #information.
        for j in range(10) :
            #index gymnastics.
            index = j + 1
            dindex = 25 + j*28

            #Calculate the surface distance between the earthquake and the
            #surface location of a vector from the center of Earth to the barycenter
            #of the given planet.
            PLat = math.radians(data[i, dindex])
            PLon = math.radians(data[i, dindex + 1])
            deltaLambda = abs(EqLon - PLon)
            temp1 = math.sin(EqLat)*math.sin(PLat)
            temp2 = math.cos(EqLat)*math.cos(PLat)*math.cos(deltaLambda)
            dist = Re*math.acos(temp1 + temp2)/1000.0 #distance in km.
            distanceEqPlanet[i, j] = dist


    #Determine the minimum distance between a earthquake and
    #a planet.
    mindistance = np.amin(distanceEqPlanet)
    print("Minimum distance {} kilometers".format(mindistance))
    
    #Now make the plot.
    planets = ['Sun','Mercury', 'Venus', 'Moon', 'Mars', 'Jupiter', 'Saturn',
               'Uranus', 'Neptune', 'Pluto']
    colors = ['skyblue', 'blue', 'green', 'red', 'brown', 'lightgray', 'yellow',
              'orange', 'hotpink', 'purple']

    #Make plots of the distributions of locations.
#    with PdfPages(outfile) as pp :
#        for i in range(0) :        
#            function_plot(distanceEqPlanet[:, i], colors[i], planets[i], pp)
    pp = 1.0
    function_plot(distanceEqPlanet[:, 3], colors[3], planets[3], pp)



############################################################
            
#Plot the force on Earth by the various planets.
def forcePlot(data) :
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages

    #Create a string holding the location of project.
    outfilepath = '/home/jdw/UM2019Autumn/M561/Project/'

    #Create a file name for the output plot.
    outfile = outfilepath + 'Plots/EqForce.pdf'
    
    #Find the number of rows in the data matrix.
    nRows = data.shape[0]

    #Pull out the force on Earth by planet during the earthquake
    #from the data array.
    dataForce = np.zeros((nRows, 10), dtype = 'float64')

    #Loop over the events and the planets to generate an array
    #of forces.
    for i in range(nRows) :
        for j in range(10) :
            index = 27 + j*28
            dataForce[i, j] = data[i, index]

    #Find the minimum value
    minForce = np.amin(dataForce)

    #Divide the dataForce array by the minimum value.
    dataForce = np.divide(dataForce, minForce)

    #Now make the plot.
    planets = ['Sun','Mercury', 'Venus', 'Moon', 'Mars', 'Jupiter', 'Saturn',
               'Uranus', 'Neptune', 'Pluto']
    colors = ['skyblue', 'blue', 'green', 'red', 'brown', 'lightgray', 'yellow',
              'orange', 'hotpink', 'purple']

    #Create a vector of evenly spaced points spanning from xmin to xmax.
    xvalues = np.linspace(1986, 2016, nRows)
    
    #Generate the figure object.
    plt.figure()
    ax = plt.axes()
    plt.grid(True)
    plt.semilogy(xvalues, dataForce[:, 0], label = planets[0])
    plt.semilogy(xvalues, dataForce[:, 1], label = planets[1])
    plt.semilogy(xvalues, dataForce[:, 2], label = planets[2])    
    plt.semilogy(xvalues, dataForce[:, 3], label = planets[3])
    plt.semilogy(xvalues, dataForce[:, 4], label = planets[4])
    plt.semilogy(xvalues, dataForce[:, 5], label = planets[5])    
    plt.semilogy(xvalues, dataForce[:, 6], label = planets[6])
    plt.semilogy(xvalues, dataForce[:, 7], label = planets[7])
    plt.semilogy(xvalues, dataForce[:, 8], label = planets[8])    
    plt.semilogy(xvalues, dataForce[:, 9], label = planets[9])

    plt.title('Plot of Normalized Force Versus Time')
    plt.xlabel('Time (years)')
    plt.ylabel('Normalized Force')
    plt.subplots_adjust(right=0.8)
    plt.legend(bbox_to_anchor = (1.04,1), borderaxespad=0)
    #Get the pdfpages object.
    pp = PdfPages(outfile)
    pp.savefig()

    #Finally, the multipage pdf object has to be closed:
    pp.close()
