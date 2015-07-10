import numpy as np

# do them all
def doall (dir = "simha/") :
    for metal in [10, 14, 17, 20, 22] :
        for start in [0.7, 1.0, 1.5, 2.0] :
            for trunc in [7, 9, 11, 13] :
                for tau in [0.3, 0.7, 1.0, 1.3, 2.0, 9.0, 13.0] :
                    for theta in [-0.175, -0.524, -0.785, -1.047, -1.396] :
                        file = "s-" + str(metal) + "-" +str(start) + "-"
                        file = file + str(trunc) + "-" + str(tau) + str(theta)
                        file1 = dir + file + "-z.mags"
                        file2 = dir + file + "-0.mags"
                        majorMinor(file1)
                        majorMinor(file2)


#   Log(Z/Zsol):  0.000
#   Fraction of blue HB stars:  0.200; Ratio of BS to HB stars:  0.000
#   Shift to TP-AGB [log(Teff),log(Lbol)]:  0.00  0.00
#   IMF: 1
#   Mag Zero Point: AB (not relevant for spec/indx files)
#   SFH: Tage= 14.96 Gyr, log(tau/Gyr)=  0.954, const=  0.000, fb=  0.000, tb=  11.00 Gyr, sf_start=  1.500, dust=(  0.00,  0.00)
#   SFH 5: tau=  9.000 Gyr, sf_start=  1.500 Gyr, sf_trunc=  9.000 Gyr, sf_theta= -0.785, dust=(  0.00,  0.00)
#
#   z log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)
#20.000  5.5000 -70.0000   0.0000 -70.0000  99.000  99.000  99.000  99.000  99.000
#20.000  5.5250 -70.0000   0.0000 -70.0000  99.000  99.000  99.000  99.000  99.000

def majorMinor ( majorFile ) :
    scale = .03   # we'll add a three percent by mass minority population

    pwords = majorFile.split("/")
    fname = pwords[-1]
    words = fname.split("-")
    minorName = "s-4"
    for j in range(2,len(words)) :
        minorName = minorName + "-" + words[j]
    minorFile = ""
    for j in range(0,len(pwords)-1) :
        minorFile = minorFile + pwords[j] + "/"
    minorFile = minorFile + minorName
    outputFile = majorFile.replace(".mags","m.mags")
    dummyNum1 = -70.0000
    dummyNum2 = 99.000
    
    majorHdr, majorData = read(majorFile)
    minorHdr, minorData = read(minorFile)
    zed, age, mass, lum, sfr, u, g, r, i, z = lineToNP(majorData) 
    mzed, mage, mmass, mlum, msfr, mu, mg, mr, mi, mz = lineToNP(minorData) 

    ix = np.nonzero((mu != dummyNum2) & (u != dummyNum2))
    u[ix] = 10**(-0.4*u[ix]) + scale*10**(-0.4*mu[ix])
    u[ix] = -2.5*np.log10(u[ix])
    ix = np.nonzero((mg != dummyNum2) & (g != dummyNum2))
    g[ix] = 10**(-0.4*g[ix]) + scale*10**(-0.4*mg[ix])
    g[ix] = -2.5*np.log10(g[ix])
    ix = np.nonzero((mr != dummyNum2) & (r != dummyNum2))
    r[ix] = 10**(-0.4*r[ix]) + scale*10**(-0.4*mr[ix])
    r[ix] = -2.5*np.log10(r[ix])
    ix = np.nonzero((mi != dummyNum2) & (i != dummyNum2))
    i[ix] = 10**(-0.4*i[ix]) + scale*10**(-0.4*mi[ix])
    i[ix] = -2.5*np.log10(i[ix])
    ix = np.nonzero((mz != dummyNum2) & (z != dummyNum2))
    z[ix] = 10**(-0.4*z[ix]) + scale*10**(-0.4*mz[ix])
    z[ix] = -2.5*np.log10(z[ix])

    ix = np.nonzero((mlum != dummyNum1) & (lum != dummyNum1))
    lum[ix] = np.log10(10**(lum[ix]) + scale*10**(mlum[ix]))
    ix = np.nonzero((msfr != dummyNum1) & (sfr != dummyNum1))
    sfr[ix] = np.log10(10**(sfr[ix]) + scale*10**(msfr[ix]))
    ix = np.nonzero((mmass != dummyNum1) & (mass != dummyNum1))
    mass[ix] = np.log10(10**(mass[ix]) + scale*10**(mmass[ix]))

    fd = open(outputFile,"w")
    fd.write(majorHdr[0])
    fd.write("#   + {:.1f}% by mass Log(Z/Zsol): {:6.3f}\n".format(
        100*scale, float(minorHdr[0][16:23])))
    for j in range(1,len(majorHdr)) :
        fd.write(majorHdr[j])
    for j in range(0,zed.size) :
        fd.write("{:6.3f} ".format( zed[j]))
        fd.write("{:7.4f} ".format(age[j]))
        fd.write("{:8.4f} ".format(mass[j]))
        fd.write("{:8.4f} ".format(lum[j]))
        fd.write("{:8.4f} ".format(sfr[j]))
        fd.write("{:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f}\n".format(
            u[j],g[j],r[j],i[j],z[j]))
    fd.close()

def read (file) :
    dummyNum = -70.0000
    header = []
    data = []
    fd = open(file, "r") 
    for line in fd :
        if line[0] == "#" :
            header.append(line)
        else :
            words = line.split()
            zed = float(words[0])
            age = float(words[1])
            mass = float(words[2])
            lum = float(words[3])
            sfr = float(words[4])
            u = float(words[5])
            g = float(words[6])
            r = float(words[7])
            i = float(words[8])
            z = float(words[9])
            if age != dummyNum :
                newLine = [zed, age, mass, lum, sfr, u, g, r, i, z]
                data.append(newLine)
    return header, data

# convert to np.arrays
def lineToNP(data) :
    zed, age, mass, lum, sfr = [],[],[],[],[]
    u, g, r, i, z = [],[],[],[],[]
    for j in range(0,len(data)) :
        zed.append(data[j][0])
        age.append(data[j][1])
        mass.append(data[j][2])
        lum.append(data[j][3])
        sfr.append(data[j][4])
        u.append(data[j][5])
        g.append(data[j][6])
        r.append(data[j][7])
        i.append(data[j][8])
        z.append(data[j][9])
    zed = np.array(zed)
    age = np.array(age)
    mass = np.array(mass)
    lum = np.array(lum)
    sfr = np.array(sfr)
    u = np.array(u)
    g = np.array(g)
    r = np.array(r)
    i = np.array(i)
    z = np.array(z)
    return zed, age, mass, lum, sfr, u, g, r, i, z
        


# end
