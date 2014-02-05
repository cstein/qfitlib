"""
Generate 
"""
vdw_radii = {}
vdw_radii[1]  = 1.2
vdw_radii[6]  = 1.5
vdw_radii[7]  = 1.5
vdw_radii[8]  = 1.4
vdw_radii[16] = 1.89

"""

"""
import numpy
EPS = 1.0e-8

def grid(ntes, r, z, rscal):
    temp_coords = []
    final_coords = []

    for i, (ri, zi) in enumerate(zip(r, z)):
        riscal = vdw_radii[zi] * rscal
        icoords = unit_sphere( ntes )
        icoords *= riscal
        for n,nc in enumerate(icoords):
            ncontact = 1
            icoords[n] = ri - icoords[n]

            for k, (rk, zk) in enumerate(zip(r, z)):
                if i == k: continue
                rkscal = vdw_radii[zk] * rscal
                dr = icoords[n] - rk
                dr2=numpy.dot( dr, dr )
                if dr2 < rkscal**2:
                    ncontact += 1

            if ncontact == 1:
                final_coords.append( icoords[n] )

    return numpy.array( final_coords )

def unit_sphere( ntes ):
    neq = int(numpy.sqrt(numpy.pi * ntes))
    nvt = neq / 2

    coordinates = []
    for i in range(nvt+1):
        angle = numpy.pi * i / nvt
        z = numpy.cos( angle )
        xy = numpy.sin( angle )
        nbo = int(xy*neq + EPS)
        if nbo < 1: nbo = 1
        for j in range(nbo):
            angle2 =  2*numpy.pi * j / nbo
            coordinates.append( numpy.array(
                  [xy*numpy.cos( angle2 ),
                   xy*numpy.sin( angle2 ),
                   z]))

    return numpy.array(coordinates)

if __name__ == '__main__':
    import sys
    import obabel
    mol = obabel.Molecule(sys.argv[1])
    R = []
    Z = []
    for atom in mol.getAtoms():
        R.append( numpy.array([atom.GetX(), atom.GetY(), atom.GetZ()]) )
        Z.append(atom.GetAtomicNum())


    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    alphas = [1.0, 0.6, 0.4, 0.2]
    colors = ['#222222','#666666','#999999','#BABABA']

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    for i in range(4):
        full_grid = grid(48, R, Z, 1.4+i*0.2)
        x,y,z = full_grid.transpose()
        ax.scatter(x,y,z,alpha=alphas[i],color=colors[i])

    plt.show()
    #plt.savefig('mep.png')
