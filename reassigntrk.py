from numpy import *
import os

def circshift(t, n):
    return hstack((t[-n:], t[:-n]))

def find(x):
    indices, = nonzero(array([x]).ravel())
    nind = indices.size
    return (indices, nind)

def sortor(x):
    ordered = x.argsort(kind='mergesort')
    return (x[ordered],ordered)

def findmax(x):
    rmax=x.argmax(0)[0]
    cmax=x.argmax(1)[0]
    return (rmax,cmax)

def unq(array):
    q = array.ravel()
    qshift = circshift(q, -1)
    indices, nind = find(q != qshift)
    if nind > 0:
        return indices
    else:
        return q.size-1
    
def mapunq(array, idx):
    #idx should be a vector array
    s = array.size
    q = array.ravel()[asarray(idx)]
    qshift = circshift(q, -1)
    indices, nind = find(q != qshift)
    if nind > 0:
        return idx[indices]
    else:
        return q.size-1


def luberize(tracks):
    #% reassigns the unique ID# to 0,1,2,3...
    #% /presort will sort on ID# first, then reassign
    #% start will begin with that ID#
    #% function returns a new track array

    ndat = tracks.shape[1]
    newtracks = tracks

    u = unq(newtracks[:,ndat-1]) +1
    ntracks = u.size

    u = r_[0,u]

    for i in r_[:ntracks]:
        newtracks[u[i]:u[i+1],ndat-1] = i  #i? or i-1?

    return newtracks

class GlobalParameters():
    istart = 1
    blocksize = 0
    ntrack = 0
    z = 0
    wp = []

    olist = zeros((1,2))    

    n = -1
    eyes = []
    pos = []
    zspan = 0
    resx = []
    bigresx = []
    mem = []
    uniqid = []
    maxid = 0
    dumphash = []
    nvalid = []
    ispan = -1
    m = -1
    xyi = []
    found = []
    si = []
    spos = []
    dimm = []
    nblocks = 1.
    map = []
    labelx = []
    labely = []

    def __init__(self):
        self.initwk()

    def initwk(self):
        self.isort = 0
        self.coltot = []
        self.rowtot = []
        self.which1 = []
        self.strt = []
        self.fnsh = []

    def cleanwk(self):
        self.initwk()

    def resetdb(self):
        self.numbonds = 0
        self.minbonds = []
        self.nclust = 0
        self.maxsz = -1
        self.mxsz = -1
        self.mysz = -1
        self.bmap = []
        self.who = []
        self.pt = []
        self.st = []
        self.fi = []
        self.hp = []
        self.ok = []
        self.unew = []
        self.nnew = []
        self.uold = []
        self.nold = []
        self.lost = []
        self.losttot = []
        self.nlost = 0
        self.lensq = 0
        self.mndisq = 0

        self.xdim = size(self.labelx)
        self.ydim = size(self.labely)
        #%  make a list of the non-trivial bonds
        self.bonds = ones((1,2))
        self.bondlen = [0]

def checktimevec(gp):
    """Check the input time vector"""
    #check the input time vector is ok, i.e. sorted and uniform

    st_t = circshift(gp.t,1)
    st = gp.t[1:size(gp.t)] - st_t[1:size(gp.t)]

    if sum(st[st<0]) != 0:
        print "The time vectors are not in order\n"
        return

    gp.wp, z = find(st>0)
    if z == 0:
        print "All positions are at the same time... go back!"
        return
    elif sum(st[gp.wp] - st[gp.wp[0]]) != 0:
        print 'WARNING - Time vector gapped or not evenly gapped!'

    gp.z = z+1

def calcblocksize(xyzs, gp):
    """calculate a blocksize which may be greater than maxdisp, but which
    keeps nblocks reasonably small."""

    volume = 1
    for d in r_[:gp.dim]:
        minn = min(xyzs[gp.wp,d])
        maxx = max(xyzs[gp.wp,d])
        volume = volume * (maxx-minn)

    blocksize = max(r_[gp.maxdisp,(volume/(20*gp.n))**(1.0/gp.dim)])
    #Tailor the factor in bottom for the particular system    

def fixzspan(gp):
    gp.zspan = 50
    if gp.n > 200:
        gp.zspan = 20
    elif gp.n > 500:
        gp.zspan = 10
    return gp.zspan

def gonext(i,xyzs,gp):
    gp.ispan = i % gp.zspan
    ##   Get the new particle positions.
    gp.m = gp.res[i+1] - gp.res[i]
    gp.eyes = gp.res[i] + r_[:gp.m] # check
    if gp.m > 0:
        gp.xyi = xyzs[gp.eyes,:gp.dim]
        gp.found = repeat(0, gp.m)

def makecube(gp):
    """construct the vertices of a 3x3x3... d-dimensional hypercube"""
    cube = zeros((3**dim, dim))
    for d in r_[:dim]:
        numb = 0
        for j in r_[0:3**dim:3**d]:
            cube[j:j+3**d,d] = numb
            numb = (numb+1) % 3
    return cube

def rastermetrictrivialbonds(gp):
    ##   construct "s", a one dimensional parameterization of the space
    ##   ( which consists of the d-dimensional raster scan of the volume.)

    abi = fix(gp.xyi/gp.blocksize)
    abpos = fix(gp.pos/gp.blocksize)
    gp.si = repeat(0, gp.m)
    gp.spos = repeat(0, gp.n)
    gp.dimm = repeat(0, gp.dim)
    gp.nblocks = 1.

    for j in r_[:gp.dim]:
        minn = min(r_[abi[:,j],abpos[:,j]])
        maxx = max(r_[abi[:,j],abpos[:,j]])
        abi[:,j] = abi[:,j] - minn
        abpos[:,j] = abpos[:,j] - minn
        gp.dimm[j] = maxx - minn + 1
        gp.si = gp.si + abi[:,j] * gp.nblocks
        gp.spos = gp.spos + abpos[:,j] * gp.nblocks
        gp.nblocks = gp.dimm[j] * gp.nblocks

def calcscoord(gp,cube):
    """ trim down (intersect) the hypercube if its too big to fit in the particle volume.
    (i.e. if dimm(j) lt 3) and then calculate the <s> coordinates of hypercube (with a corner @
    the origin) shift the hypercube <s> coordinates to be centered around the origin"""

    cub = cube;
    deg, ndeg = find(gp.dimm<3)
    if ndeg != 0:
        for j in r_[:deg.size]:
            cub = cub[nonzero(cub[:,deg[j]] < gp.dimm[deg[j]]), :]

    scube = repeat(0, cube.shape[0])
    coff = 1
    for j in r_[:gp.dim]:
        scube = scube + cube[:,j] * coff
        coff = coff * gp.dimm[j]

    coff = 1
    for j in r_[:gp.dim]:
        if gp.dimm[j] > 3:
            scube = scube - coff
        coff = gp.dimm[j] * coff

    scube = (scube + gp.nblocks) % gp.nblocks
    return scube

def prepmat(gp):
    mg = mgrid[:gp.ntrack,:gp.m]
    xmat = mg[1]
    ymat = mg[0]
    return xmat, ymat

def update(gp, xyzs):
    gp.wp, nww = find(gp.resx[gp.ispan,:] >= 0)
    if nww > 0:
        gp.pos[gp.wp,:] = xyzs[gp.resx[gp.ispan,gp.wp].astype('int64'), :gp.dim]
        if gp.goodenough > 0:
            gp.nvalid[gp.wp] = gp.nvalid[gp.wp] + 1
    else:
        print 'Warning, tracking zero particles!\n'

    #we need to add new guys, as appropriate.
    newguys, nnew = find(gp.found == 0)

    if nnew > 0:
        newarr = zeros((gp.zspan, nnew)) - 1
        gp.resx = c_[gp.resx, newarr]
        gp.resx[gp.ispan, gp.n:] = gp.eyes[newguys]
        gp.pos = vstack([gp.pos, xyzs[gp.eyes[newguys],:gp.dim]])
        gp.mem = r_[gp.mem, repeat(0, nnew)]
        gp.uniqid = r_[gp.uniqid, r_[:nnew] + gp.maxid]
        gp.maxid = gp.maxid + nnew
        if gp.goodenough > 0:
            gp.dumphash = r_[gp.dumphash, repeat(0,nnew)]
            gp.nvalid = r_[gp.nvalid, repeat(1,nnew)]

        gp.n = gp.n + nnew

def updatemem(gp, i):
    ##   update the 'memory' array
    gp.wp, nok = find(gp.resx[gp.ispan,:] != -1)

    if nok != 0:
        gp.mem[gp.wp] = 0
    gp.mem = gp.mem + (gp.resx[gp.ispan,:].T == -1)

    wlost, nlost = find(gp.mem == gp.memory + 1)
    if nlost > 0:
        gp.pos[wlost,:] = - gp.maxdisp
        if gp.goodenough > 0:
            wdump, ndump = find(gp.nvalid[wlost] < gp.goodenough)
            if ndump > 0:
                gp.dumphash[wlost[wdump]] = 1

    if (gp.ispan == gp.zspan-1) or (i== gp.z-1):
        nnew = gp.n - gp.bigresx.shape[1]
        if nnew > 0:
            newarr = zeros((gp.z, nnew)) -1
            gp.bigresx = column_stack([gp.bigresx, newarr])

        if gp.goodenough > 0:
            if sum(gp.dumphash) > 0:
                wkeep, nkeep = find(gp.dumphash == 0)
                gp.resx = gp.resx[:,wkeep]
                gp.bigresx = gp.bigresx[:,wkeep]
                gp.pos = gp.pos[wkeep,:]
                gp.mem = gp.mem[wkeep]
                gp.uniqid = gp.uniqid[wkeep]
                gp.nvalid = gp.nvalid[wkeep]
                gp.n = nkeep
                gp.dumphash = repeat(0, nkeep)

        if gp.quiet != 1:
            print "{:d} of {:d} done. Tracking {:d} particles, {:d} tracks total\n".format(i, gp.z, ntrk, gp.n)

        gp.bigresx[i-gp.ispan:i+1,:] = gp.resx[:gp.ispan+1,:]
        gp.resx = zeros((gp.zspan,gp.n)) -1

        wpull, npull = find(gp.pos[:,0] == -gp.maxdisp)

        if npull > 0:
            lillist = zeros((1,2))
            for ipull in r_[:npull]:
                wpull2, npull2 = find(gp.bigresx[:,wpull[ipull]] != -1)
                thing = column_stack((gp.bigresx[wpull2, wpull[ipull]], repeat(0, npull2) + gp.uniqid[wpull[ipull]]))
                lillist = vstack([lillist, thing])

            gp.olist = vstack((gp.olist,lillist[1:,:]))

        wkeep, nkeep = find(gp.pos[:,0] >= 0)

        if nkeep == 0:
            print 'Were going to crash now, no particles....\n'
        gp.resx = gp.resx[:,wkeep]
        gp.bigresx = gp.bigresx[:,wkeep]
        gp.pos = gp.pos[wkeep,:]
        gp.mem = gp.mem[wkeep]
        gp.uniqid = gp.uniqid[wkeep]
        gp.n = nkeep
        gp.dumphash = repeat(0, nkeep)
        if gp.goodenough > 0:
            gp.nvalid = gp.nvalid[wkeep]

def initfornotnsqrd(gp):
    ed, gp.isort = sortor(gp.si)

    gp.strt, gp.fnsh = gp.settargetrange()

    gp.coltot = repeat(0, gp.m)
    gp.rowtot = repeat(0, gp.n)
    gp.which1 = repeat(0, gp.n)

def settargetrange(gp):
    """make a hash table which will allow us to know which new particles
    are at a given si."""

    gp.strt = repeat(0, gp.nblocks) -1    # -1 is a tag
    gp.fnsh = repeat(0, gp.nblocks)

    for j in r_[:gp.m]:
        if gp.strt[gp.si[gp.isort[j]]] == -1:    # if it is the first time for filling the block
            gp.strt[gp.si[gp.isort[j]]] = j      # the beginning of gp.si with that value.
            gp.fnsh[gp.si[gp.isort[j]]] = j
        else:
            gp.fnsh[gp.si[gp.isort[j]]] = j      # only the finish is being updated while gp.si takes the same value.
    return (gp.strt, gp.fnsh)

def findtrivialbonds (gp, j):
    """ find those trivial bonds"""
    distq = repeat(0, gp.map.size)
    for d in r_[:gp.dim]:
        distq = distq + (gp.xyi[gp.map,d] - gp.pos[j,d])**2
    return distq

def countgood(gp):
    gp.rowtot = repeat(0, gp.n)
    gp.rowtot[gp.wh] = sum(gp.ltmax,axis=1)
    if gp.ntrack > 1:
        gp.coltot = repeat(0, gp.ntrack)            
        gp.coltot = sum(gp.ltmax,axis=0)
    else:
        gp.coltot = gp.ltmax

def labelxy(gp):
    gp.wp, ngood = find(gp.rowtot == 1)
    if ngood != 0:
        ww, ngood = find(gp.coltot[gp.which1[gp.wp]] == 1)
        if ngood != 0:
            gp.resx[gp.ispan,gp.wp[ww]] = gp.eyes[gp.which1[gp.wp[ww]]]
            gp.found[gp.which1[gp.wp[ww]]] = 1
            gp.rowtot[gp.wp[ww]] = 0
            gp.coltot[gp.which1[gp.wp[ww]]] = 0

    gp.labely, ngood = find(gp.rowtot>0)
    if ngood != 0:
        gp.labelx = find(gp.coltot > 0)[0]
        nontrivial = True
    else:
        nontrivial = False

    return nontrivial

def indexmax(gp):
    gp.which1 = repeat(0, gp.n)
    for j in r_[:gp.ntrack]:
        gp.wp = gp.ltmax[j,:].argmax()
        if gp.wp.size > 1: gp.wp = gp.wp[0]
        gp.which1[gp.wh[j]] = gp.wp

    return gp.which1

def checkpt(gp):
    if gp.pt[gp.who] != gp.st[gp.who] - 1:            
        gp.ok[gp.hp[gp.pt[gp.who]]] = 1            

def extractsubnetwork(gp):
    """Extracts connected sub-networks of the non-trivial bonds.
    NB: lista/b can have redundant entries due to multiple-connected subnetworks"""

    lista = array(repeat(1,gp.numbonds))
    listb = array(repeat(1,gp.numbonds))

    gp.nclust = 0
    gp.maxsz = 0
    thru = gp.xdim

    while thru != 0:
        gp.wp = find(gp.bonds[:,1] >= 0)[0]
        lista[0] = array(gp.bonds[gp.wp[0],1])
        listb[0] = array(gp.bonds[gp.wp[0],0])
        gp.bonds[gp.wp[0],:] = - (gp.nclust+1)      # consider -1 as a tag "Done"
        adda = 1
        addb = 1
        donea = 0    # donea and doneb serve as indices later
        doneb = 0
        #         adda  = 2; addb  = 2;
        #         donea = 1; doneb = 1;

        cycle = True
        while (cycle):
            if donea != adda:                    
                gp.wp, ngood = find(gp.bonds[:,1] == lista[donea])
                if ngood != 0:
                    listb[addb:addb+ngood] = gp.bonds[gp.wp,0]
                    gp.bonds[gp.wp,:] = - (gp.nclust+1)     # consider -1 as a tag "Done"
                    addb = addb + ngood
                donea = donea + 1
            if doneb != addb:                    
                gp.wp, ngood = find(gp.bonds[:,0] == listb[doneb])
                if ngood != 0:
                    lista[adda:adda+ngood] = gp.bonds[gp.wp,1]
                    gp.bonds[gp.wp,:] = - (gp.nclust+1)       # consider -1 as a tag "Done"
                    adda = adda + ngood
                doneb = doneb + 1

            if donea == adda and doneb == addb:                    
                cycle = False
            else:
                cycle = True

        tempb, idxsortb = sortor(listb[:doneb])
        tempa, idxsorta = sortor(lista[:donea])
        xsz = size( mapunq(listb[:doneb], idxsortb) )
        ysz = size( mapunq(lista[:donea], idxsorta) )

        if xsz*ysz > gp.maxsz:
            gp.maxsz = xsz*ysz
            gp.mxsz = xsz
            gp.mysz = ysz

        thru = thru - xsz
        gp.nclust = gp.nclust + 1

    gp.bmap = gp.bonds[:,0]

    return

def updatebonds(gp):
    if gp.who == gp.nnew-1:
        ww = find(gp.lost == 0)[0]
        dsq = sum(gp.lensq[gp.pt[ww]]) + gp.losttot * gp.maxdisq
        if dsq < gp.mndisq:
            gp.minbonds = gp.pt[ww]
            gp.mndisq = dsq
    else:
        gp.who = gp.who +1

def evallost(gp):
    notlost = -gp.lost[gp.who] -1
    if (notlost % 2 == 1) and (gp.losttot != gp.nlost):
        gp.lost[gp.who] = 1
        gp.losttot = gp.losttot + 1
        checkpt(gp)
        updatebonds(gp)
    else:
        checkpt(gp)
        gp.pt[gp.who] = gp.st[gp.who] - 1
        if gp.lost[gp.who]:
            gp.lost[gp.who] = 0
            gp.losttot = gp.losttot - 1
        gp.who = gp.who - 1


def trackmem(xyzs, maxdisp=5, memory=3, dim=2, goodenough=3, quiet=1):

    xyzs=xyzs[xyzs[:,-1].argsort(kind='merge'),:]
    gp = GlobalParameters()

    gp.maxdisp = maxdisp
    gp.memory = memory
    gp.dim = dim
    gp.goodenough = goodenough
    gp.maxdisq = maxdisp ** 2
    gp.dd = xyzs.shape[1]
    gp.quiet = quiet
    gp.t = xyzs[:,-1]
    checktimevec(gp)
    gp.res = r_[0, unq(gp.t)+1,gp.t.size]

    #res indexes the border of time frames
    gp.n = gp.res[1] - gp.res[0]
    gp.eyes = r_[:gp.n]
    gp.pos = xyzs[gp.eyes,:gp.dim]
    #Cut out x,y position data spanninng first period
    gp.zspan = fixzspan(gp)
    gp.resx = zeros((gp.zspan, gp.n)) -1
    gp.bigresx = zeros((gp.z, gp.n)) -1
    gp.mem = repeat(0,gp.n)
    gp.uniqid = r_[:gp.n]
    gp.maxid = gp.n
    gp.resx[0,:] = gp.eyes
    # eyes = 0,1,2,...(first break point)
    ##resx(1,:) = eyes;     #put the first set of feature indices in the first row of resx
    if gp.goodenough > 0:
        #goodenough is a parameter regarding length
        gp.dumphash = repeat(0, gp.n)
        gp.nvalid = repeat(1, gp.n)

    notnsqrd = 1*(gp.n > 200 and gp.dim < 7)
    ##Use fancy code for large n, small d
    ##notnsqrd = (sqrt(n*ngood) > 200) & (dim < 7);

    if notnsqrd:
        cube = makecube()
        calcblocksize(xyzs, gp)

    for i in r_[gp.istart:gp.z]:
        gonext(i,xyzs,gp)        
        if gp.m > 0:
        #   THE TRIVIAL BOND CODE BEGINS
            if notnsqrd:
                rastermetrictrivialbonds(gp)
                scube = calcscoord(cube,gp)
                gp.initwk()
                initfornotnsqrd(gp)
                for j in r_[:gp.n]:
                    gp.map = -1
                    s = (scube + gp.spos[j]) % gp.nblocks
                    gp.wp, ngood = find(gp.strt[s] != -1)
                    if ngood != 0:
                        s = s[gp.wp]
                        for k in r_[:ngood]:
                            gp.map = r_[gp.map, gp.isort[gp.strt[s[k]]:gp.fnsh[s[k]]]]
                        gp.map = gp.map[1:]
                        distq = findtrivialbonds(gp, j)
                        gp.wp, gp.rowtot[j] = find(distq < gp.maxdisq)
                        if gp.rowtot[j] > 0:
                            gp.coltot[gp.map[gp.wp]] = gp.coltot[gp.map[gp.wp]] + 1
                            gp.which1[j] = gp.map[gp.wp[0]]

                nontrivial = labelxy(gp)
                gp.cleanwk()
                ##clear abi,clear abpos,clear fnsh, clear rowtot, clear coltot, clear which1, clear isort

            else:
                gp.initwk()
                gp.wh, gp.ntrack = find(gp.pos[:,0] >= 0)
                if gp.ntrack == 0:
                    print "There are no valid particles to track idiot!"
                    break
                xmat, ymat = prepmat(gp)

                for d in r_[:gp.dim]:
                    x = gp.xyi[:,d]                    
                    y = gp.pos[gp.wh,d]

                    if d == 0:
                        dq = (x[xmat] - y[ymat.T].T)**2
                    else:
                        dq = dq + (x[xmat] - y[ymat.T].T)**2
                    
                gp.ltmax = (dq < gp.maxdisq)
                #% figure out which trivial bonds go with which        
                countgood(gp)
                gp.which1 = indexmax(gp)
                ntrk = fix(gp.n - sum(gp.rowtot==0))

                nontrivial = labelxy(gp)

                gp.cleanwk()
                #            del rowtot, coltot, which1

            if nontrivial:
                gp.resetdb()
                for j in r_[:gp.ydim]:
                    distq = repeat(0, gp.xdim)
                    for d in r_[:gp.dim]:
                        distq = distq + (gp.xyi[gp.labelx,d].ravel() - gp.pos[gp.labely[j],d].ravel())**2
                    gp.wp, ngood = find(distq < gp.maxdisq)
                    gp.bonds = vstack((gp.bonds,column_stack((gp.wp+1,repeat(1, ngood)+j))))
                    gp.bondlen = r_[gp.bondlen,distq[gp.wp]]

                gp.bonds = gp.bonds[1:,:]
                gp.bondlen = gp.bondlen.ravel()[1:]
                gp.numbonds = gp.bonds[:,0].size
                mbonds = gp.bonds.copy()       ### This is the pitfall of python !!!!!
            
                if max(r_[gp.xdim,gp.ydim]) < 4:    #xdim and ydim are related to size.
                    gp.nclust = 1
                    gp.maxsz = 0
                    gp.mxsz = gp.xdim
                    gp.mysz = gp.ydim
                    gp.bmap = repeat(0, gp.bonds.shape[0]) -1
                else:
                    #%   THE SUBNETWORK CODE BEGINS
                    extractsubnetwork(gp)
                    #% THE SUBNETWORK CODE ENDS

                #%   THE PERMUTATION CODE BEGINS

                for nc in r_[:gp.nclust]:
                    gp.wp, gp.nbonds = find(gp.bmap == -1*(nc+1))

                    gp.bonds = mbonds[gp.wp,:]
                    gp.lensq = gp.bondlen[gp.wp]
                    temp, sortidx = sortor(gp.bonds[:,0])
                    gp.uold = array([gp.bonds[ mapunq(gp.bonds[:,0],sortidx), 0 ]]).ravel()
                    gp.nold = gp.uold.size
                    gp.unew = array(gp.bonds[ unq(gp.bonds[:,1]), 1 ]).ravel()
                    gp.nnew = gp.unew.size
                    
                    # check that runtime is not excessive
                    if gp.nnew > 5:
                        rnsteps = 1
                        for ii in r_[:gp.nnew]:
                            rnsteps = rnsteps * find(gp.bonds[:,1] == gp.unew[ii])[1]
                            if rnsteps > 5.e+4:
                                print ' Warning: difficult combinatorics encountered.\n'
                                print ' Program may not finish- Try reducing maxdisp.\n'
                            if rnsteps > 2.e+5:
                                print ' Excessive Combinatorics! Try reducing maxdisp.\n'
                                return

                    gp.st = repeat(0, gp.nnew)
                    gp.fi = repeat(0, gp.nnew)
                    gp.hp = repeat(0, gp.nbonds)
                    gp.ok = repeat(1, gp.nold) +1
                    if gp.nnew-gp.nold > 0:
                        gp.nlost = gp.nnew - gp.nold
                    else:
                        gp.nlost=0;

                    for ii in r_[:gp.nold]:
                        gp.hp[find(gp.bonds[:,0] == gp.uold[ii])[0]] = ii
                    gp.st[0] = 0
                    gp.fi[gp.nnew-1] = gp.nbonds -1
                    if gp.nnew > 1:
                        sb = gp.bonds[:,1]
                        sbr = circshift(sb,1)
                        sbl = circshift(sb,-1)
                        gp.st[1:] = find(sb[1:] != sbr[1:] )[0] + 1
                        gp.fi[:gp.nnew-1] = find(sb[:gp.nbonds-1] != sbl[:gp.nbonds-1])[0]

                    checkflag = 0
                    while checkflag != 2:
                        gp.pt = gp.st -1
                        gp.lost = repeat(0, gp.nnew)
                        gp.who = 0
                        gp.losttot = 0
                        gp.mndisq = gp.nnew * gp.maxdisq
                        while gp.who != -1:
                            if gp.pt[gp.who] != gp.fi[gp.who]:
                                gp.wp, ngood = find(gp.ok[gp.hp[gp.pt[gp.who]+1:gp.fi[gp.who]+1]])
                                if ngood > 0:
                                    checkpt(gp)
                                    gp.pt[gp.who] = gp.pt[gp.who] + gp.wp[0] +1
                                    gp.ok[gp.hp[gp.pt[gp.who]]] = 0
                                    updatebonds(gp)
                                else:
                                    evallost(gp)                                    
                            else:
                                evallost(gp)

                        checkflag = checkflag + 1
                        if checkflag == 1:
                        #   we need to check that our constraint on nlost is not forcing us away from the minimum id's
                            plost = min(r_[fix([gp.mndisq/gp.maxdisq]), gp.nnew-1])
                            if plost > gp.nlost + 1:
                                gp.nlost = plost
                            else:
                                checkflag = 2
                                
                    #%   update resx using the minimum bond configuration
                    gp.resx[gp.ispan, gp.labely[gp.bonds[gp.minbonds,1].astype('int64')-1]] = gp.eyes[gp.labelx[gp.bonds[gp.minbonds,0].astype('int64')-1]]
                    gp.found[gp.labelx[gp.bonds[gp.minbonds,0].astype('int64')-1]] = 1

                    #                del db

                #%   THE PERMUTATION CODE ENDS
            #     here we want to update our initial position estimates
            update(gp,xyzs)
        else:
            print ' Warning- No positions found for t='

        updatemem(gp, i)

    ## the big loop over z time steps....

    #%%  make a final scan for short trajectories that weren't lost at the end.
    if gp.goodenough > 0:
        nvalid = (gp.bigresx > 0).sum(axis=0)
        wkeep, nkeep = find( nvalid >= gp.goodenough)
        if nkeep < gp.n:
            gp.bigresx = gp.bigresx[:,wkeep]
            gp.n = nkeep;
            gp.uniqid = gp.uniqid[wkeep];
            gp.pos = gp.pos[wkeep,:];

    #%  make the final scan to 'pull' everybody else into the olist.

    wpull, npull = find( gp.pos[:,0] != -2*gp.maxdisp )
    if npull > 0:

        lillist = array([1,1])
        for ipull in r_[:npull]:
            wpull2, npull2 = find(gp.bigresx[:,wpull[ipull]] != -1)
            lillist = vstack([ lillist, column_stack([ gp.bigresx[wpull2,wpull[ipull]],repeat(0,npull2) + gp.uniqid[wpull[ipull]] ] )])

        gp.olist = vstack([gp.olist, lillist[1:,:]])
        
    gp.olist = gp.olist[1:,:]
    
#  free up a little memory for the final step!
    gp.bigresx = 0
    gp.resx = 0

    gp.res = zeros((gp.olist.shape[0],gp.dd+1))

    for j in r_[:gp.dd]:
        gp.res[:,j] = xyzs[gp.olist[:,0].astype('int64'),j]

    gp.res[:,gp.dd] = gp.olist[:,1]
    
    if gp.res.shape[0] > 0:
        lub = luberize(gp.res)
    else:
        lub = gp.res

    # % end of uberize code
    return lub

def concatenatepolished(folderpath="PolishedSpots"):
    """Concatenate Polished Spots files. Experimental"""
    xyzs = zeros([1,6])
    #    xyzs = zeros([1,3])    
    for file in os.listdir(folderpath):
        if file.startswith('PSpA'):
            frame=int(file[4:8])
            filepath=folderpath+'/'+file
            data=genfromtxt(filepath,skiprows=1)
            data=data[data[:,12]==1,:]
            #            data=data[data[:,12]==1,4:6]
            #            xyzs=row_stack([xyzs,column_stack([data,repeat(frame,data.shape[0])])])
    xyzs=row_stack([xyzs,column_stack([data[:,1:3],data[:,4:6],data[:,0],data[3]])])
    savetxt("xyzs.txt",xyzs[1:,],fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%d')
    return xyzs[1:,:]

def concatenatespotfiles(folderpath="OriginalSpots"):
    """Concatenate Original Spots files"""    
    xyzs = zeros([1,3])
    for file in os.listdir(folderpath):
        frame=int(file[3:7])
        filepath=folderpath+'/'+file
        data=genfromtxt(filepath,skiprows=1)
        xyzs=row_stack([xyzs,column_stack([data[:,1:3],repeat(frame,data.shape[0])])])
    savetxt("xyzs.txt",xyzs[1:,],fmt='%10.5f\t%10.5f\t%d')        
    return xyzs[1:,:]

def gapfind(t):
    return find(asarray([i not in t for i in r_[:max(t)+1]]))

def interpolategaps(data):
    idmax=data[data.shape[0]-1,data.shape[1]-1]
    print idmax
    olist = repeat(0.0,4)    
    for cid in r_[:idmax+1]:
        #        print "cid=",cid
        pind,=where(data[:,data.shape[1]-1]==cid)
        subdata=data[pind,:]
        subind=subdata[:,2] - subdata[0,2]
        pdata=zeros((max(subind)+1,4))
        for i in r_[:max(subind)+1]:
            if i in subind:
                pdata[i,:]=subdata[subind==i,:]
            else:
                headind=nonzero([j < i for j in subind])
                headlast=max(subind[headind])
                tailind,=nonzero([j > i for j in subind])
                tailfirst=min(subind[tailind])
                gaplength=tailfirst - headlast
                alpha=(i-headlast)*1.0/(tailfirst-headlast)
                pdata[i,:] = r_[subdata[subind==headlast,0]*(1.0-alpha)+subdata[subind==tailfirst,0]*alpha, subdata[subind==headlast,1]*(1.0-alpha)+subdata[subind==tailfirst,1]*alpha,i+subdata[0,2],cid]
        olist = vstack([olist,pdata])
      
    return olist[1:,:]

def rearrangecolumns(data):
    return column_stack([data[:,3],data[:,2],data[:,:2]])



def main():
#    xyzs = concatenatespotfiles('PolishedSpots')
    xyzs = concatenatespotfiles()
    newtracks=trackmem(xyzs, dim=2, memory=3)
    #    newtracks=trackmem(xyzs,maxdisp=5, memory=3, dim=2, goodenough=3)    
    savetxt("opl_newtrk.txt",rearrangecolumns(newtracks),fmt='%d\t%d\t%10.5f\t%10.5f')
    interpolated=interpolategaps(newtracks)
    newtracks=rearrangecolumns(interpolated)
    savetxt("opl_reass.txt",newtracks, fmt='%d\t%d\t%10.5f\t%10.5f')

if __name__ == '__main__':
    main()
