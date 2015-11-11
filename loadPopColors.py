import numpy as np

# do them all
#   lib can be miles or basel
def doall ( dir = "simha/", lib="basel") :
    data = []
    metallicites_basel = [10,14,17, 20, 22]
    metallicites_miles = [2,3,4, 5]        
    metallicites = metallicites_basel
    if lib == "miles" : metallicites = metallicites_miles
    for metal in metallicites:
        #for start in [0.7, 1.0, 1.5, 2.0] :
        #    for trunc in [7, 9, 11, 13] :
        #        for tau in [0.3, 1.0, 1.3, 2.0, 9.0, 13.0] :
        for start in [1.5, 2.0] :
            for trunc in [11, 13] :
                for tau in [1.3, 2.0, 9.0, 13.0] :
        #for start in [0.7, 1.0 ] :
        #    for trunc in [7, 9] :
        #        for tau in [0.3, 1.0] :
                    for theta in [-0.175, -0.524, -0.785, -1.047, -1.396] :
                        file = "s-" + str(metal) + "-" +str(start) + "-"
                        file = file + str(trunc) + "-" + str(tau) + str(theta)
                        file = dir + file + ".mags"
                        ug,gr,ri,iz,grr,gir,kii,kri,ml = loadSplines(file)
                        data.append([ug,gr,ri,iz,grr,gir,kii,kri,ml])
                        
    return data

#   lib can be miles or basel
def doallnames ( dir = "simha/", lib="basel") :
    data = []
    metallicites_basel = [10,14,17, 20, 22]
    metallicites_miles = [2,3,4,5]        
    metallicites = metallicites_basel
    if lib == "miles" : metallicites = metallicites_miles
    counter = 0
    names = dict()
    for metal in metallicites:
        for start in [0.7, 1.0, 1.5, 2.0] :
            for trunc in [7, 9, 11, 13] :
                for tau in [0.3, 1.0, 1.3, 2.0, 9.0, 13.0] :
                    for theta in [-0.175, -0.524, -0.785, -1.047, -1.396] :
                        file = "s-" + str(metal) + "-" +str(start) + "-"
                        file = file + str(trunc) + "-" + str(tau) + str(theta)
                        file = dir + file + ".mags"
                        names[counter] = file       
                        names[counter,"metals"] = metal
                        counter += 1
    return names


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

def loadSplines ( filename ) :
    from scipy.interpolate import interp1d
    from scipy.interpolate import UnivariateSpline

    dummyNum1 = -70.0000
    dummyNum2 = 99.000
    
    hdr, data = read(filename)
    zed, age, mass, lum, sfr, ug, gr, ri, iz, grr, gir, kii, kri = lineToNP(data) 

    ix = np.nonzero( (zed > 0) & (zed < 1.6))
    ug = UnivariateSpline(zed[ix], ug[ix], s=0)
    gr = UnivariateSpline(zed[ix], gr[ix], s=0)
    ri = UnivariateSpline(zed[ix], ri[ix], s=0)
    iz = UnivariateSpline(zed[ix], iz[ix], s=0)
    grr = UnivariateSpline(zed[ix], grr[ix], s=0)
    gir = UnivariateSpline(zed[ix], gir[ix], s=0)
    kii = UnivariateSpline(zed[ix], kii[ix], s=0)
    kri = UnivariateSpline(zed[ix], kri[ix], s=0)
    ml = UnivariateSpline(zed[ix], mass[ix]-lum[ix], s=0)

    return ug,gr,ri,iz,grr,gir,kii,kri,ml
    

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
            ug = float(words[5])
            gr = float(words[6])
            ri = float(words[7])
            iz = float(words[8])
            grr = float(words[9])
            gir = float(words[10])
            kii = float(words[11])
            kri = float(words[12])
            if age != dummyNum :
                newLine = [zed, age, mass, lum, sfr, ug, gr, ri, iz, grr,gir,kii,kri]
                data.append(newLine)
    return header, data

# convert to np.arrays
def lineToNP(data) :
    zed, age, mass, lum, sfr = [],[],[],[],[]
    ug, gr, ri, iz, grr, gir = [],[],[],[],[],[]
    kii, kri = [],[]
    for j in range(len(data)-1,-1,-1) :
        zed.append(data[j][0])
        age.append(data[j][1])
        mass.append(data[j][2])
        lum.append(data[j][3])
        sfr.append(data[j][4])
        ug.append(data[j][5])
        gr.append(data[j][6])
        ri.append(data[j][7])
        iz.append(data[j][8])
        grr.append(data[j][9])
        gir.append(data[j][10])
        kii.append(data[j][11])
        kri.append(data[j][12])
    zed = np.array(zed)
    age = np.array(age)
    mass = np.array(mass)
    lum = np.array(lum)
    sfr = np.array(sfr)
    ug = np.array(ug)
    gr = np.array(gr)
    ri = np.array(ri)
    iz = np.array(iz)
    grr = np.array(grr)
    gir = np.array(gir)
    kii = np.array(kii) ;# k-correction taking obs i to restframe i: kii= i_o - i_obs
    kri = np.array(kri) ;# k-correction taking obs i to restframe r: kri= r_o - i_obs
    return zed, age, mass, lum, sfr, ug, gr, ri, iz, grr, gir, kii, kri
        
# end
