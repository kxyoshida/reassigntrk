from numpy import *
import os

def concatenatespotfiles(folderpath):
    """Concatenate Polished Spots files"""
    xyzs = zeros([1,3])
    for file in os.listdir(folderpath):
        if file.startswith('PSpA'):
            frame=int(file[4:8])
            filepath=folderpath+'/'+file
            data=genfromtxt(filepath,skiprows=1)
            data=data[data[:,12]==1,4:6]
            xyzs=row_stack([xyzs,column_stack([data,repeat(frame,data.shape[0])])])
    savetxt("xyzs.txt",xyzs[1:,],fmt='%10.5f\t%10.5f\t%d')
    return xyzs[1:,:]

class parameters:
    """A class defining global parameters"""
    dim = 2
    maxdisp = 5
    maxdisq = maxdisp ** 2
    goodenough  =   3    #param.good
    memory_b = 3
    #    This has set to 3 for AP2GFP
    #    memory_b = 1    
    quiet = 1
    istart = 1
    blocksize = 0
    ntrack = 0
    dd = 0
    t = []
    z = 0
    res = []
    wp = []
    olist = []

    def __init__(self, xyzs):
        self.dd = xyzs.shape[1]
        self.t = xyzs[:,-1]
        self.checktimevec()
        self.res = r_[0, unq(self.t)+1,self.t.size]
        # unq always returns the border indices - 1 so that the value should be
        # compensated by adding 1. Also the last term could be the last indice of self.t,
        # and 1<=unq<=self.t.size-1, i.e. 2<=unq+1<=self.t.size.
        # Actually the last two terms become the same and I do not know how it could be justified.
        ##res = unq(t,[]);
        ##res = res+1;
        ##res = [1,res,length(t)+1];

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


    def calcblocksize(self, xyzs,n):
        """calculate a blocksize which may be greater than maxdisp, but which
        keeps nblocks reasonably small."""

        volume = 1
        for d in r_[:dim]:
            minn = min(xyzs[self.wp,d])
            maxx = max(xyzs[self.wp,d])
            volume = volume * (maxx-minn)

        blocksize = max(r_[maxdisp,(volume/(20*n))**(1.0/dim)])
        #Tailor the factor in bottom for the particular system    
        

class segment:
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

    def __init__(self, xyzs, par):
        #res indexes the border of time frames
        self.n = par.res[1] - par.res[0]
        self.eyes = r_[:self.n]
        self.pos = xyzs[self.eyes,:par.dim]
        #Cut out x,y position data spanninng first period
        self.zspan = self.fixzspan()
        self.resx = zeros((self.zspan, self.n)) -1
        self.bigresx = zeros((par.z, self.n)) -1
        self.mem = repeat(0,self.n)
        self.uniqid = r_[:self.n]
        self.maxid = self.n
        self.resx[0,:] = self.eyes
        # eyes = 0,1,2,...(first break point)
        ##resx(1,:) = eyes;     #put the first set of feature indices in the first row of resx
        if par.goodenough > 0:
            #goodenough is a parameter regarding length
            self.dumphash = repeat(0, self.n)
            self.nvalid = repeat(1, self.n)

    def fixzspan(self):
        self.zspan = 50
        if self.n > 200:
            self.zspan = 20
        elif self.n > 500:
            self.zspan = 10
        return self.zspan

    def gonext(self,i,par,xyzs):
        self.ispan = i % self.zspan
        ##   Get the new particle positions.
        self.m = par.res[i+1] - par.res[i]
        self.eyes = par.res[i] + r_[:self.m] # check
        if self.m > 0:
            self.xyi = xyzs[self.eyes,:par.dim]
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

    def rastermetrictrivialbonds(self, par):
        ##   construct "s", a one dimensional parameterization of the space
        ##   ( which consists of the d-dimensional raster scan of the volume.)

        abi = fix(self.xyi/par.blocksize)
        abpos = fix(self.pos/par.blocksize)
        self.si = repeat(0, self.m)
        self.spos = repeat(0, self.n)
        self.dimm = repeat(0, par.dim)
        self.nblocks = 1.

        for j in r_[:par.dim]:
            minn = min(r_[abi[:,j],abpos[:,j]])
            maxx = max(r_[abi[:,j],abpos[:,j]])
            abi[:,j] = abi[:,j] - minn
            abpos[:,j] = abpos[:,j] - minn
            self.dimm[j] = maxx - minn + 1
            self.si = self.si + abi[:,j] * self.nblocks
            self.spos = self.spos + abpos[:,j] * self.nblocks
            self.nblocks = self.dimm[j] * self.nblocks

    def calcscoord(self,cube,par):
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
        for j in r_[:par.dim]:
            scube = scube + cube[:,j] * coff
            coff = coff * self.dimm[j]

        coff = 1
        for j in r_[:par.dim]:
            if self.dimm[j] > 3:
                scube = scube - coff
            coff = self.dimm[j] * coff

        scube = (scube + self.nblocks) % self.nblocks
        return scube

    def prepmat(self, par):
        mg = mgrid[:par.ntrack,:self.m]
        xmat = mg[1]
        ymat = mg[0]
        return xmat, ymat

    def update(self, par, xyzs):
        par.wp, nww = find(self.resx[self.ispan,:] >= 0)
        if nww > 0:
            self.pos[par.wp,:] = xyzs[self.resx[self.ispan,par.wp].astype('int64'), :par.dim]
            if par.goodenough > 0:
                self.nvalid[par.wp] = self.nvalid[par.wp] + 1
        else:
            print 'Warning, tracking zero particles!\n'

        #we need to add new guys, as appropriate.
        newguys, nnew = find(self.found == 0)

        if nnew > 0:
            newarr = zeros((self.zspan, nnew)) - 1
            self.resx = c_[self.resx, newarr]
            self.resx[self.ispan, self.n:] = self.eyes[newguys]
            self.pos = vstack([self.pos, xyzs[self.eyes[newguys],:par.dim]])
            self.mem = r_[self.mem, repeat(0, nnew)]
            self.uniqid = r_[self.uniqid, r_[:nnew] + self.maxid]
            self.maxid = self.maxid + nnew
            if par.goodenough > 0:
                self.dumphash = r_[self.dumphash, repeat(0,nnew)]
                self.nvalid = r_[self.nvalid, repeat(1,nnew)]

            self.n = self.n + nnew

    def updatemem(self, par, i):
        ##   update the 'memory' array
        par.wp, nok = find(self.resx[self.ispan,:] != -1)
        
        if nok != 0:
            self.mem[par.wp] = 0
        self.mem = self.mem + (self.resx[self.ispan,:].T == -1)

        wlost, nlost = find(self.mem == par.memory_b + 1)
        if nlost > 0:
            self.pos[wlost,:] = - par.maxdisp
            if par.goodenough > 0:
                wdump, ndump = find(self.nvalid[wlost] < par.goodenough)
                if ndump > 0:
                    self.dumphash[wlost[wdump]] = 1

        if (self.ispan == self.zspan-1) or (i== par.z-1):
            nnew = self.n - self.bigresx.shape[1]
            if nnew > 0:
                newarr = zeros((par.z, nnew)) -1
                self.bigresx = column_stack([self.bigresx, newarr])

            if par.goodenough > 0:
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

            if par.quiet != 1:
                print "{:d} of {:d} done. Tracking {:d} particles, {:d} tracks total\n".format(i, par.z, ntrk, self.n)

            self.bigresx[i-self.ispan:i+1,:] = self.resx[:self.ispan+1,:]
            self.resx = zeros((self.zspan,self.n)) -1

            wpull, npull = find(self.pos[:,0] == -par.maxdisp)
            
            if npull > 0:
                lillist = zeros((1,2))
                for ipull in r_[:npull]:
                    wpull2, npull2 = find(self.bigresx[:,wpull[ipull]] != -1)
                    thing = column_stack((self.bigresx[wpull2, wpull[ipull]], repeat(0, npull2) + self.uniqid[wpull[ipull]]))
                    lillist = vstack([lillist, thing])

                par.olist = vstack((par.olist,lillist[1:,:]))

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
            if par.goodenough > 0:
                self.nvalid = self.nvalid[wkeep]

class work:
    hp = []
    isort = 0
    coltot = []
    rowtot = []
    which1 = []
    strt = []
    fnsh = []

    def initfornotnsqrd(self, seg):
        ed, self.isort = sortor(seg.si)

        self.hp = nonzero(seg.si == 0)
        self.strt, self.fnsh = self.settargetrange(seg)

        self.coltot = repeat(0, seg.m)
        self.rowtot = repeat(0, seg.n)
        self.which1 = repeat(0, seg.n)

    def settargetrange(self, seg):
        """make a hash table which will allow us to know which new particles
        are at a given si."""

        self.strt = repeat(0, seg.nblocks) -1    # -1 is a tag
        self.fnsh = repeat(0, seg.nblocks)
    
        for j in r_[:seg.m]:
            if self.strt[seg.si[self.isort[j]]] == -1:    # if it is the first time for filling the block
                self.strt[seg.si[self.isort[j]]] = j      # the beginning of seg.si with that value.
                self.fnsh[seg.si[self.isort[j]]] = j
            else:
                self.fnsh[seg.si[self.isort[j]]] = j      # only the finish is being updated while seg.si takes the same value.
        return (self.strt, self.fnsh)

    def findtrivialbonds (self, seg, par, j):
        """ find those trivial bonds"""
        distq = repeat(0, self.map.size)
        for d in r_[:par.dim]:
            distq = distq + (seg.xyi[self.map,d] - seg.pos[j,d])**2
        return distq

    def countgood(self, ltmax, wh, par, seg):
        self.rowtot = repeat(0, seg.n)
        self.rowtot[wh] = sum(ltmax,axis=1)
        if par.ntrack > 1:
            self.coltot = repeat(0, par.ntrack)            
            self.coltot = sum(ltmax,axis=0)
        else:
            self.coltot = ltmax
            
    def labelxy(self, seg, par):
        par.wp, ngood = find(self.rowtot == 1)
        if ngood != 0:
            ww, ngood = find(self.coltot[self.which1[par.wp]] == 1)
            if ngood != 0:
                seg.resx[seg.ispan,par.wp[ww]] = seg.eyes[self.which1[par.wp[ww]]]
                seg.found[self.which1[par.wp[ww]]] = 1
                self.rowtot[par.wp[ww]] = 0
                self.coltot[self.which1[par.wp[ww]]] = 0
            
        seg.labely, ngood = find(self.rowtot>0)
        if ngood != 0:
            seg.labelx = find(self.coltot > 0)[0]
            nontrivial = True
        else:
            nontrivial = False
        
        return nontrivial

    def indexmax(self,ltmax,wh,par,seg):
        self.which1 = repeat(0, seg.n)
        for j in r_[:par.ntrack]:
            par.wp = ltmax[j,:].argmax()
            if par.wp.size > 1: par.wp = par.wp[0]
            self.which1[wh[j]] = par.wp

        return self.which1
    

class database:
    xdim = 0
    ydim = 0
    bonds = []
    bondlen = []
    numbonds = 0
    minbonds = []
    nclust = 0
    maxsz = -1
    mxsz = -1
    mysz = -1
    bmap = []
    who = []
    pt = []
    st = []
    fi = []
    hp = []
    ok = []
    unew = []
    nnew = []
    uold = []
    nold = []
    lost = []
    losttot = []
    nlost = 0
    lensq = 0
    mndisq = 0

    def __init__(self, seg):
        self.xdim = size(seg.labelx)
        self.ydim = size(seg.labely)
        #%  make a list of the non-trivial bonds
        self.bonds = ones((1,2))
        self.bondlen = [0]

    def checkpt(self):
        if self.pt[self.who] != self.st[self.who] - 1:            
            self.ok[self.hp[self.pt[self.who]]] = 1            

    def extractsubnetwork(self, par):
        """Extracts connected sub-networks of the non-trivial bonds.
        NB: lista/b can have redundant entries due to multiple-connected subnetworks"""

        lista = array(repeat(1,self.numbonds))
        listb = array(repeat(1,self.numbonds))
    
        self.nclust = 0
        self.maxsz = 0
        thru = self.xdim

        while thru != 0:
            par.wp = find(self.bonds[:,1] >= 0)[0]
            lista[0] = array(self.bonds[par.wp[0],1])
            listb[0] = array(self.bonds[par.wp[0],0])
            self.bonds[par.wp[0],:] = - (self.nclust+1)      # consider -1 as a tag "Done"
            adda = 1
            addb = 1
            donea = 0    # donea and doneb serve as indices later
            doneb = 0
            #         adda  = 2; addb  = 2;
            #         donea = 1; doneb = 1;

            cycle = True
            while (cycle):
                if donea != adda:                    
                    par.wp, ngood = find(self.bonds[:,1] == lista[donea])
                    if ngood != 0:
                        listb[addb:addb+ngood] = self.bonds[par.wp,0]
                        self.bonds[par.wp,:] = - (self.nclust+1)     # consider -1 as a tag "Done"
                        addb = addb + ngood
                    donea = donea + 1
                if doneb != addb:                    
                    par.wp, ngood = find(self.bonds[:,0] == listb[doneb])
                    if ngood != 0:
                        lista[adda:adda+ngood] = self.bonds[par.wp,1]
                        self.bonds[par.wp,:] = - (self.nclust+1)       # consider -1 as a tag "Done"
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


    def updatebonds(self, par):
        if self.who == self.nnew-1:
            ww = find(self.lost == 0)[0]
            dsq = sum(self.lensq[self.pt[ww]]) + self.losttot * par.maxdisq
            if dsq < self.mndisq:
                self.minbonds = self.pt[ww]
                self.mndisq = dsq
        else:
            self.who = self.who +1

    def evallost(self, par):
        notlost = -self.lost[self.who] -1
        if (notlost % 2 == 1) and (self.losttot != self.nlost):
            self.lost[self.who] = 1
            self.losttot = self.losttot + 1
            self.checkpt()
            self.updatebonds(par)
        else:
            self.checkpt()
            self.pt[self.who] = self.st[self.who] - 1
            if self.lost[self.who]:
                self.lost[self.who] = 0
                self.losttot = self.losttot - 1
            self.who = self.who - 1

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

    
def track(xyzs):

    xyzs=xyzs[xyzs[:,-1].argsort(kind='merge'),:]

    par = parameters(xyzs)
    seg = segment(xyzs,par)
    par.olist = zeros((1,2))
    ##olist = [1,1];

    notnsqrd = 1*(seg.n > 200 and par.dim < 7)
    ##Use fancy code for large n, small d
    ##notnsqrd = (sqrt(n*ngood) > 200) & (dim < 7);

    if notnsqrd:
        cube = par.makecube()
        par.calcblocksize(xyzs, par.wp, seg.n)

    for i in r_[par.istart:par.z]:
        seg.gonext(i,par,xyzs)
        if seg.m > 0:
        #   THE TRIVIAL BOND CODE BEGINS
            if notnsqrd:
                seg.rastermetrictrivialbonds(par)
                scube = calcscoord(cube,par)
                wk = work()
                wk.initfornotnsqrd(seg)
                for j in r_[:seg.n]:
                    seg.map = -1
                    s = (scube + seg.spos[j]) % seg.nblocks
                    par.wp, ngood = find(wk.strt[s] != -1)
                    if ngood != 0:
                        s = s[par.wp]
                        for k in r_[:ngood]:
                            seg.map = r_[seg.map, wk.isort[wk.strt[s[k]]:wk.fnsh[s[k]]]]
                        seg.map = seg.map[1:]
                        distq = seg.findtrivialbonds(seg,par, j)
                        par.wp, wk.rowtot[j] = find(distq < par.maxdisq)
                        if wk.rowtot[j] > 0:
                            wk.coltot[seg.map[par.wp]] = wk.coltot[seg.map[par.wp]] + 1
                            wk.which1[j] = seg.map[par.wp[0]]

                nontrivial = wk.labelxy(seg, par)
                del wk
                ##clear abi,clear abpos,clear fnsh, clear rowtot, clear coltot, clear which1, clear isort

            else:
                wk = work()                
                wh, par.ntrack = find(seg.pos[:,0] >= 0)
                if par.ntrack == 0:
                    print "There are no valid particles to track idiot!"
                    break
                xmat, ymat = seg.prepmat(par)

                for d in r_[:par.dim]:
                    x = seg.xyi[:,d]                    
                    y = seg.pos[wh,d]

                    if d == 0:
                        dq = (x[xmat] - y[ymat.T].T)**2
                    else:
                        dq = dq + (x[xmat] - y[ymat.T].T)**2
                    
                ltmax = (dq < par.maxdisq)
                #% figure out which trivial bonds go with which        
                wk.countgood(ltmax,wh,par,seg)
                wk.which1 = wk.indexmax(ltmax,wh,par,seg)
                ntrk = fix(seg.n - sum(wk.rowtot==0))

                nontrivial = wk.labelxy(seg, par)

                del wk
                #            del rowtot, coltot, which1

            if nontrivial:
                db = database(seg)
                for j in r_[:db.ydim]:
                    distq = repeat(0, db.xdim)
                    for d in r_[:par.dim]:
                        distq = distq + (seg.xyi[seg.labelx,d].ravel() - seg.pos[seg.labely[j],d].ravel())**2
                    par.wp, ngood = find(distq < par.maxdisq)
                    db.bonds = vstack((db.bonds,column_stack((par.wp+1,repeat(1, ngood)+j))))
                    db.bondlen = r_[db.bondlen,distq[par.wp]]

                db.bonds = db.bonds[1:,:]
                db.bondlen = db.bondlen.ravel()[1:]
                db.numbonds = db.bonds[:,0].size
                mbonds = db.bonds.copy()       ### This is the pitfall of python !!!!!
            
                if max(r_[db.xdim,db.ydim]) < 4:    #xdim and ydim are related to size.
                    db.nclust = 1
                    db.maxsz = 0
                    db.mxsz = db.xdim
                    db.mysz = db.ydim
                    db.bmap = repeat(0, db.bonds.shape[0]) -1
                else:
                    #%   THE SUBNETWORK CODE BEGINS
                    db.extractsubnetwork(par)
                    #% THE SUBNETWORK CODE ENDS

                #%   THE PERMUTATION CODE BEGINS

                for nc in r_[:db.nclust]:
                    par.wp, db.nbonds = find(db.bmap == -1*(nc+1))

                    db.bonds = mbonds[par.wp,:]
                    db.lensq = db.bondlen[par.wp]
                    temp, sortidx = sortor(db.bonds[:,0])
                    db.uold = array([db.bonds[ mapunq(db.bonds[:,0],sortidx), 0 ]]).ravel()
                    db.nold = db.uold.size
                    db.unew = array(db.bonds[ unq(db.bonds[:,1]), 1 ]).ravel()
                    db.nnew = db.unew.size
                    
                    # check that runtime is not excessive
                    if db.nnew > 5:
                        rnsteps = 1
                        for ii in r_[:db.nnew]:
                            rnsteps = rnsteps * find(db.bonds[:,1] == db.unew[ii])[1]
                            if rnsteps > 5.e+4:
                                print ' Warning: difficult combinatorics encountered.\n'
                                print ' Program may not finish- Try reducing maxdisp.\n'
                            if rnsteps > 2.e+5:
                                print ' Excessive Combinatorics! Try reducing maxdisp.\n'
                                return

                    db.st = repeat(0, db.nnew)
                    db.fi = repeat(0, db.nnew)
                    db.hp = repeat(0, db.nbonds)
                    db.ok = repeat(1, db.nold) +1
                    if db.nnew-db.nold > 0:
                        db.nlost = db.nnew - db.nold
                    else:
                        db.nlost=0;

                    for ii in r_[:db.nold]:
                        db.hp[find(db.bonds[:,0] == db.uold[ii])[0]] = ii
                    db.st[0] = 0
                    db.fi[db.nnew-1] = db.nbonds -1
                    if db.nnew > 1:
                        sb = db.bonds[:,1]
                        sbr = circshift(sb,1)
                        sbl = circshift(sb,-1)
                        db.st[1:] = find(sb[1:] != sbr[1:] )[0] + 1
                        db.fi[:db.nnew-1] = find(sb[:db.nbonds-1] != sbl[:db.nbonds-1])[0]

                    checkflag = 0
                    while checkflag != 2:
                        db.pt = db.st -1
                        db.lost = repeat(0, db.nnew)
                        db.who = 0
                        db.losttot = 0
                        db.mndisq = db.nnew * par.maxdisq
                        while db.who != -1:
                            if db.pt[db.who] != db.fi[db.who]:
                                par.wp, ngood = find(db.ok[db.hp[db.pt[db.who]+1:db.fi[db.who]+1]])
                                if ngood > 0:
                                    db.checkpt()
                                    db.pt[db.who] = db.pt[db.who] + par.wp[0] +1
                                    db.ok[db.hp[db.pt[db.who]]] = 0
                                    db.updatebonds(par)
                                else:
                                    db.evallost(par)
                            else:
                                db.evallost(par)

                        checkflag = checkflag + 1
                        if checkflag == 1:
                        #   we need to check that our constraint on nlost is not forcing us away from the minimum id's
                            plost = min(r_[fix([db.mndisq/par.maxdisq]), db.nnew-1])
                            if plost > db.nlost + 1:
                                db.nlost = plost
                            else:
                                checkflag = 2
                                
                    #%   update resx using the minimum bond configuration
                    seg.resx[seg.ispan, seg.labely[db.bonds[db.minbonds,1].astype('int64')-1]] = seg.eyes[seg.labelx[db.bonds[db.minbonds,0].astype('int64')-1]]
                    seg.found[seg.labelx[db.bonds[db.minbonds,0].astype('int64')-1]] = 1

                    #                del db

                #%   THE PERMUTATION CODE ENDS
            #     here we want to update our initial position estimates
            seg.update(par,xyzs)
        else:
            print ' Warning- No positions found for t='

        seg.updatemem(par, i)

    ## the big loop over z time steps....

    #%%  make a final scan for short trajectories that weren't lost at the end.
    if par.goodenough > 0:
        nvalid = (seg.bigresx > 0).sum(axis=0)
        wkeep, nkeep = find( nvalid >= par.goodenough)
        if nkeep < seg.n:
            seg.bigresx = seg.bigresx[:,wkeep]
            seg.n = nkeep;
            seg.uniqid = seg.uniqid[wkeep];
            seg.pos = seg.pos[wkeep,:];

    #%  make the final scan to 'pull' everybody else into the olist.

    wpull, npull = find( seg.pos[:,0] != -2*par.maxdisp )
    if npull > 0:

        lillist = array([1,1])
        for ipull in r_[:npull]:
            wpull2, npull2 = find(seg.bigresx[:,wpull[ipull]] != -1)
            lillist = vstack([ lillist, column_stack([ seg.bigresx[wpull2,wpull[ipull]],repeat(0,npull2) + seg.uniqid[wpull[ipull]] ] )])

        par.olist = vstack([par.olist, lillist[1:,:]])
        
    par.olist = par.olist[1:,:]
    
#  free up a little memory for the final step!
    seg.bigresx = 0
    seg.resx = 0

    par.res = zeros((par.olist.shape[0],par.dd+1))

    for j in r_[:par.dd]:
        par.res[:,j] = xyzs[par.olist[:,0].astype('int64'),j]

    par.res[:,par.dd] = par.olist[:,1]
    
    if par.res.shape[0] > 0:
        lub = luberize(par.res)
    else:
        lub = res

    # % end of uberize code
    return lub

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
    xyzs = concatenatespotfiles('PolishedSpots')
    newtracks=track(xyzs)
    savetxt("opl_newtrk.txt",rearrangecolumns(newtracks),fmt='%d\t%d\t%10.5f\t%10.5f')
    interpolated=interpolategaps(newtracks)
    newtracks=rearrangecolumns(interpolated)
    savetxt("opl_reass.txt",newtracks, fmt='%d\t%d\t%10.5f\t%10.5f')

if __name__ == '__main__':
    main()
