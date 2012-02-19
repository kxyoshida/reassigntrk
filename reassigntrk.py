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

    def checktimevec(self):
        """Check the input time vector"""
        #check the input time vector is ok, i.e. sorted and uniform

        st_t = circshift(self.t,1)
        st = self.t[1:size(self.t)] - st_t[1:size(self.t)]

        if sum(st[st<0]) != 0:
            print "The time vectors are not in order\n"
            return

        self.wp, z = find(st>0)
        if z == 0:
            print "All positions are at the same time... go back!"
            return
        elif sum(st[self.wp] - st[self.wp[0]]) != 0:
            print 'WARNING - Time vector gapped or not evenly gapped!'

        self.z = z+1

    def calcblocksize(self, xyzs):
        """calculate a blocksize which may be greater than maxdisp, but which
        keeps nblocks reasonably small."""

        volume = 1
        for d in r_[:self.dim]:
            minn = min(xyzs[self.wp,d])
            maxx = max(xyzs[self.wp,d])
            volume = volume * (maxx-minn)

        blocksize = max(r_[self.maxdisp,(volume/(20*self.n))**(1.0/self.dim)])
        #Tailor the factor in bottom for the particular system    

    def fixzspan(self):
        self.zspan = 50
        if self.n > 200:
            self.zspan = 20
        elif self.n > 500:
            self.zspan = 10
        return self.zspan

    def gonext(self,i,xyzs):
        self.ispan = i % self.zspan
        ##   Get the new particle positions.
        self.m = self.res[i+1] - self.res[i]
        self.eyes = self.res[i] + r_[:self.m] # check
        if self.m > 0:
            self.xyi = xyzs[self.eyes,:self.dim]
            self.found = repeat(0, self.m)

    def makecube(self):
        """construct the vertices of a 3x3x3... d-dimensional hypercube"""
        cube = zeros((3**dim, dim))
        for d in r_[:dim]:
            numb = 0
            for j in r_[0:3**dim:3**d]:
                cube[j:j+3**d,d] = numb
                numb = (numb+1) % 3
        return cube

    def rastermetrictrivialbonds(self):
        ##   construct "s", a one dimensional parameterization of the space
        ##   ( which consists of the d-dimensional raster scan of the volume.)

        abi = fix(self.xyi/self.blocksize)
        abpos = fix(self.pos/self.blocksize)
        self.si = repeat(0, self.m)
        self.spos = repeat(0, self.n)
        self.dimm = repeat(0, self.dim)
        self.nblocks = 1.

        for j in r_[:self.dim]:
            minn = min(r_[abi[:,j],abpos[:,j]])
            maxx = max(r_[abi[:,j],abpos[:,j]])
            abi[:,j] = abi[:,j] - minn
            abpos[:,j] = abpos[:,j] - minn
            self.dimm[j] = maxx - minn + 1
            self.si = self.si + abi[:,j] * self.nblocks
            self.spos = self.spos + abpos[:,j] * self.nblocks
            self.nblocks = self.dimm[j] * self.nblocks

    def calcscoord(self,cube):
        """ trim down (intersect) the hypercube if its too big to fit in the particle volume.
        (i.e. if dimm(j) lt 3) and then calculate the <s> coordinates of hypercube (with a corner @
        the origin) shift the hypercube <s> coordinates to be centered around the origin"""

        cub = cube;
        deg, ndeg = find(self.dimm<3)
        if ndeg != 0:
            for j in r_[:deg.size]:
                cub = cub[nonzero(cub[:,deg[j]] < self.dimm[deg[j]]), :]

        scube = repeat(0, cube.shape[0])
        coff = 1
        for j in r_[:self.dim]:
            scube = scube + cube[:,j] * coff
            coff = coff * self.dimm[j]

        coff = 1
        for j in r_[:self.dim]:
            if self.dimm[j] > 3:
                scube = scube - coff
            coff = self.dimm[j] * coff

        scube = (scube + self.nblocks) % self.nblocks
        return scube

    def prepmat(self):
        mg = mgrid[:self.ntrack,:self.m]
        xmat = mg[1]
        ymat = mg[0]
        return xmat, ymat

    def update(self, xyzs):
        self.wp, nww = find(self.resx[self.ispan,:] >= 0)
        if nww > 0:
            self.pos[self.wp,:] = xyzs[self.resx[self.ispan,self.wp].astype('int64'), :self.dim]
            if self.goodenough > 0:
                self.nvalid[self.wp] = self.nvalid[self.wp] + 1
        else:
            print 'Warning, tracking zero particles!\n'

        #we need to add new guys, as appropriate.
        newguys, nnew = find(self.found == 0)

        if nnew > 0:
            newarr = zeros((self.zspan, nnew)) - 1
            self.resx = c_[self.resx, newarr]
            self.resx[self.ispan, self.n:] = self.eyes[newguys]
            self.pos = vstack([self.pos, xyzs[self.eyes[newguys],:self.dim]])
            self.mem = r_[self.mem, repeat(0, nnew)]
            self.uniqid = r_[self.uniqid, r_[:nnew] + self.maxid]
            self.maxid = self.maxid + nnew
            if self.goodenough > 0:
                self.dumphash = r_[self.dumphash, repeat(0,nnew)]
                self.nvalid = r_[self.nvalid, repeat(1,nnew)]

            self.n = self.n + nnew

    def updatemem(self, i):
        ##   update the 'memory' array
        self.wp, nok = find(self.resx[self.ispan,:] != -1)

        if nok != 0:
            self.mem[self.wp] = 0
        self.mem = self.mem + (self.resx[self.ispan,:].T == -1)

        wlost, nlost = find(self.mem == self.memory + 1)
        if nlost > 0:
            self.pos[wlost,:] = - self.maxdisp
            if self.goodenough > 0:
                wdump, ndump = find(self.nvalid[wlost] < self.goodenough)
                if ndump > 0:
                    self.dumphash[wlost[wdump]] = 1

        if (self.ispan == self.zspan-1) or (i== self.z-1):
            nnew = self.n - self.bigresx.shape[1]
            if nnew > 0:
                newarr = zeros((self.z, nnew)) -1
                self.bigresx = column_stack([self.bigresx, newarr])

            if self.goodenough > 0:
                if sum(self.dumphash) > 0:
                    wkeep, nkeep = find(self.dumphash == 0)
                    self.resx = self.resx[:,wkeep]
                    self.bigresx = self.bigresx[:,wkeep]
                    self.pos = self.pos[wkeep,:]
                    self.mem = self.mem[wkeep]
                    self.uniqid = self.uniqid[wkeep]
                    self.nvalid = self.nvalid[wkeep]
                    self.n = nkeep
                    self.dumphash = repeat(0, nkeep)

            if self.quiet != 1:
                print "{:d} of {:d} done. Tracking {:d} particles, {:d} tracks total\n".format(i, self.z, ntrk, self.n)

            self.bigresx[i-self.ispan:i+1,:] = self.resx[:self.ispan+1,:]
            self.resx = zeros((self.zspan,self.n)) -1

            wpull, npull = find(self.pos[:,0] == -self.maxdisp)

            if npull > 0:
                lillist = zeros((1,2))
                for ipull in r_[:npull]:
                    wpull2, npull2 = find(self.bigresx[:,wpull[ipull]] != -1)
                    thing = column_stack((self.bigresx[wpull2, wpull[ipull]], repeat(0, npull2) + self.uniqid[wpull[ipull]]))
                    lillist = vstack([lillist, thing])

                self.olist = vstack((self.olist,lillist[1:,:]))

            wkeep, nkeep = find(self.pos[:,0] >= 0)

            if nkeep == 0:
                print 'Were going to crash now, no particles....\n'
            self.resx = self.resx[:,wkeep]
            self.bigresx = self.bigresx[:,wkeep]
            self.pos = self.pos[wkeep,:]
            self.mem = self.mem[wkeep]
            self.uniqid = self.uniqid[wkeep]
            self.n = nkeep
            self.dumphash = repeat(0, nkeep)
            if self.goodenough > 0:
                self.nvalid = self.nvalid[wkeep]

    def initfornotnsqrd(self):
        ed, self.isort = sortor(self.si)

        self.strt, self.fnsh = self.settargetrange()

        self.coltot = repeat(0, self.m)
        self.rowtot = repeat(0, self.n)
        self.which1 = repeat(0, self.n)

    def settargetrange(self):
        """make a hash table which will allow us to know which new particles
        are at a given si."""

        self.strt = repeat(0, self.nblocks) -1    # -1 is a tag
        self.fnsh = repeat(0, self.nblocks)

        for j in r_[:self.m]:
            if self.strt[self.si[self.isort[j]]] == -1:    # if it is the first time for filling the block
                self.strt[self.si[self.isort[j]]] = j      # the beginning of self.si with that value.
                self.fnsh[self.si[self.isort[j]]] = j
            else:
                self.fnsh[self.si[self.isort[j]]] = j      # only the finish is being updated while self.si takes the same value.
        return (self.strt, self.fnsh)

    def findtrivialbonds (self, j):
        """ find those trivial bonds"""
        distq = repeat(0, self.map.size)
        for d in r_[:self.dim]:
            distq = distq + (self.xyi[self.map,d] - self.pos[j,d])**2
        return distq

    def countgood(self):
        self.rowtot = repeat(0, self.n)
        self.rowtot[self.wh] = sum(self.ltmax,axis=1)
        if self.ntrack > 1:
            self.coltot = repeat(0, self.ntrack)            
            self.coltot = sum(self.ltmax,axis=0)
        else:
            self.coltot = self.ltmax

    def labelxy(self):
        self.wp, ngood = find(self.rowtot == 1)
        if ngood != 0:
            ww, ngood = find(self.coltot[self.which1[self.wp]] == 1)
            if ngood != 0:
                self.resx[self.ispan,self.wp[ww]] = self.eyes[self.which1[self.wp[ww]]]
                self.found[self.which1[self.wp[ww]]] = 1
                self.rowtot[self.wp[ww]] = 0
                self.coltot[self.which1[self.wp[ww]]] = 0

        self.labely, ngood = find(self.rowtot>0)
        if ngood != 0:
            self.labelx = find(self.coltot > 0)[0]
            nontrivial = True
        else:
            nontrivial = False

        return nontrivial

    def indexmax(self):
        self.which1 = repeat(0, self.n)
        for j in r_[:self.ntrack]:
            self.wp = self.ltmax[j,:].argmax()
            if self.wp.size > 1: self.wp = self.wp[0]
            self.which1[self.wh[j]] = self.wp

        return self.which1

    def checkpt(self):
        if self.pt[self.who] != self.st[self.who] - 1:            
            self.ok[self.hp[self.pt[self.who]]] = 1            

    def extractsubnetwork(self):
        """Extracts connected sub-networks of the non-trivial bonds.
        NB: lista/b can have redundant entries due to multiple-connected subnetworks"""

        lista = array(repeat(1,self.numbonds))
        listb = array(repeat(1,self.numbonds))

        self.nclust = 0
        self.maxsz = 0
        thru = self.xdim

        while thru != 0:
            self.wp = find(self.bonds[:,1] >= 0)[0]
            lista[0] = array(self.bonds[self.wp[0],1])
            listb[0] = array(self.bonds[self.wp[0],0])
            self.bonds[self.wp[0],:] = - (self.nclust+1)      # consider -1 as a tag "Done"
            adda = 1
            addb = 1
            donea = 0    # donea and doneb serve as indices later
            doneb = 0

            cycle = True
            while (cycle):
                if donea != adda:                    
                    self.wp, ngood = find(self.bonds[:,1] == lista[donea])
                    if ngood != 0:
                        listb[addb:addb+ngood] = self.bonds[self.wp,0]
                        self.bonds[self.wp,:] = - (self.nclust+1)     # consider -1 as a tag "Done"
                        addb = addb + ngood
                    donea = donea + 1
                if doneb != addb:                    
                    self.wp, ngood = find(self.bonds[:,0] == listb[doneb])
                    if ngood != 0:
                        lista[adda:adda+ngood] = self.bonds[self.wp,1]
                        self.bonds[self.wp,:] = - (self.nclust+1)       # consider -1 as a tag "Done"
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

            if xsz*ysz > self.maxsz:
                self.maxsz = xsz*ysz
                self.mxsz = xsz
                self.mysz = ysz

            thru = thru - xsz
            self.nclust = self.nclust + 1

        self.bmap = self.bonds[:,0]

        return

    def updatebonds(self):
        if self.who == self.nnew-1:
            ww = find(self.lost == 0)[0]
            dsq = sum(self.lensq[self.pt[ww]]) + self.losttot * self.maxdisq
            if dsq < self.mndisq:
                self.minbonds = self.pt[ww]
                self.mndisq = dsq
        else:
            self.who = self.who +1

    def evallost(self):
        notlost = -self.lost[self.who] -1
        if (notlost % 2 == 1) and (self.losttot != self.nlost):
            self.lost[self.who] = 1
            self.losttot = self.losttot + 1
            self.checkpt()
            self.updatebonds()
        else:
            self.checkpt()
            self.pt[self.who] = self.st[self.who] - 1
            if self.lost[self.who]:
                self.lost[self.who] = 0
                self.losttot = self.losttot - 1
            self.who = self.who - 1


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
    gp.checktimevec()
    gp.res = r_[0, unq(gp.t)+1,gp.t.size]

    #res indexes the border of time frames
    gp.n = gp.res[1] - gp.res[0]
    gp.eyes = r_[:gp.n]
    gp.pos = xyzs[gp.eyes,:gp.dim]
    #Cut out x,y position data spanninng first period
    gp.zspan = gp.fixzspan()
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
        gp.calcblocksize(xyzs)

    for i in r_[gp.istart:gp.z]:
        gp.gonext(i,xyzs)        
        if gp.m > 0:
        #   THE TRIVIAL BOND CODE BEGINS
            if notnsqrd:
                gp.rastermetrictrivialbonds()
                scube = gp.calcscoord(cube)
                gp.initwk()
                gp.initfornotnsqrd()
                for j in r_[:gp.n]:
                    gp.map = -1
                    s = (scube + gp.spos[j]) % gp.nblocks
                    gp.wp, ngood = find(gp.strt[s] != -1)
                    if ngood != 0:
                        s = s[gp.wp]
                        for k in r_[:ngood]:
                            gp.map = r_[gp.map, gp.isort[gp.strt[s[k]]:gp.fnsh[s[k]]]]
                        gp.map = gp.map[1:]
                        distq = gp.findtrivialbonds(j)
                        gp.wp, gp.rowtot[j] = find(distq < gp.maxdisq)
                        if gp.rowtot[j] > 0:
                            gp.coltot[gp.map[gp.wp]] = gp.coltot[gp.map[gp.wp]] + 1
                            gp.which1[j] = gp.map[gp.wp[0]]

                nontrivial = gp.labelxy()
                gp.cleanwk()
                ##clear abi,clear abpos,clear fnsh, clear rowtot, clear coltot, clear which1, clear isort

            else:
                gp.initwk()
                gp.wh, gp.ntrack = find(gp.pos[:,0] >= 0)
                if gp.ntrack == 0:
                    print "There are no valid particles to track idiot!"
                    break
                xmat, ymat = gp.prepmat()

                for d in r_[:gp.dim]:
                    x = gp.xyi[:,d]                    
                    y = gp.pos[gp.wh,d]

                    if d == 0:
                        dq = (x[xmat] - y[ymat.T].T)**2
                    else:
                        dq = dq + (x[xmat] - y[ymat.T].T)**2
                    
                gp.ltmax = (dq < gp.maxdisq)
                #% figure out which trivial bonds go with which        
                gp.countgood()
                gp.which1 = gp.indexmax()
                ntrk = fix(gp.n - sum(gp.rowtot==0))

                nontrivial = gp.labelxy()

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
                    gp.extractsubnetwork()
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
                                    gp.checkpt()
                                    gp.pt[gp.who] = gp.pt[gp.who] + gp.wp[0] +1
                                    gp.ok[gp.hp[gp.pt[gp.who]]] = 0
                                    gp.updatebonds()
                                else:
                                    gp.evallost()                                    
                            else:
                                gp.evallost()

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
            gp.update(xyzs)
        else:
            print ' Warning- No positions found for t='

        gp.updatemem(i)

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
