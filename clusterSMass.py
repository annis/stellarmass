import numpy as np
def haloStellarMass ( filename="vt-stellar-masses.txt", outfile="vt-halos-stellar-masses.txt", verbose=1 ) :

    ra,dec,z,zerr,id,central,hostid, gr,gre, gi,gie, kri,krie, kii,kiie, \
        i,distmod,rabs,iabs, mcMass, taMass, mass, std = np.genfromtxt(filename, unpack=True)

    id = np.array(id).astype(int)
    central = np.array(central).astype(int)
    hostid = np.array(hostid).astype(int)

    fd = open(outfile,"w")

    halos = np.unique(hostid)
    halos = np.sort(halos)
    if verbose: print "# hostid  median(ra) median(dec) median(z) ngals log(stellar_mass)  (h=0.7, Om=0.3, flat)"
    fd.write("# hostid  median(ra) median(dec) median(z) ngals log(stellar_mass)  (h=0.7, Om=0.3, flat)\n")
    for halo in halos:
        ix = np.nonzero(hostid == halo)
        zmedian = np.median(z[ix])
        ramedian = np.median(ra[ix])
        decmedian = np.median(dec[ix])
        ngals = ra[ix].size
        linear_mass = 10**mass[ix]
        mass_errors = np.log(10.)*linear_mass*std[ix]
        mass_std = np.sqrt((mass_errors**2).sum())

        sum_mass = linear_mass.sum()
        sum_mass_std = mass_std/sum_mass/np.log(10.)
        sum_mass = np.log10(sum_mass)
        if verbose: 
            print "{:10d} {:10.6f} {:10.5f}  {:6.3f}    {:4d}      {:6.3f} {:6.4f}".format(
            halo,ramedian, decmedian, zmedian, ngals, sum_mass, sum_mass_std)
        fd.write("{:10d} {:10.6f} {:10.5f}  {:6.3f}     {:4d}      {:6.3f} {:6.4f}\n".format(
            halo,ramedian, decmedian, zmedian, ngals, sum_mass, sum_mass_std))
    fd.close()
# cat vt-halos-stellar-masses.txt | sort -k6 -nr > vtt
