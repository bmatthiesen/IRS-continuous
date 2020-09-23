import numpy as np
import scipy.constants
import scipy.optimize
import itertools as it
import matplotlib.pyplot as plt
import pandas as pd
#from mpl_toolkits.mplot3d import Axes3D

# constants
c0 = scipy.constants.speed_of_light
er = 6371e3 # earth radius [m]
GM =  3.986004418e14 # geocentric gravitational constant ("mu" in document) [m^3 s^−2]

# Rotation matrix around x axis
def Rx(a):
    return np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])

def _reshape(p):
    if p.ndim == 1:
        if p.shape[0] != 3:
            raise ValueError("not an Euclidean vector")

        return p[:,np.newaxis]
    elif p.ndim ==2 :
        if p.shape[1] == 3:
            return p.T
        elif p.shape[0] == 3:
            return p
        else:
            raise ValueError("not an Euclidean vector")
    else:
        raise ValueError("only 1D and 2D arrays supported")


# rot: rotation of IRS plane
# dtx: horizontal distance of Tx to IRS (i.e., in z direction) [m]
# drx: horizontal distance of Tx to IRS (i.e., in z direction) [m]
# elevation of IRS over ground (Tx is on ground) [m]
# satH: satellite height above ground
#fc: Frequency [Hz]
# Tpoints: density of time interval
# visAng: minimum Tx-Rx angle for transmission [deg]
class SatIRS:
    def __init__(self, M, N, dx, dy = None, rot = 0, dtx = 1000, drx = 1000, fc = 2e9, h = 100, satH = 1500e3, mu = 1, isoIRS = False, Tpoints = 100, visAng = 10, T = None):
        if dy is None:
            dy = dx

        # set effective size of iso IRS
        if isoIRS:
            Airs = M*N*dx*dy # physical size of IRS
            Aiso = (c0 / fc)**2 / 4 / np.pi # effective area of isotropic antenna
            aspect = M/N

            self.N = np.sqrt(Airs / Aiso / aspect)
            self.M = int(np.round(aspect * self.N))
            self.N = int(np.round(self.N))

            self.dx = M*dx / self.M
            self.dy = N*dy / self.N
        else:
            self.M = M
            self.N = N
            self.dx = dx
            self.dy = dy


        # save parameter
        self.rotation = rot
        self.fc = fc
        self.mu = mu
        self.isoIRS = isoIRS

        # Tx position
        self.pt = Rx(rot) @ np.array([0, -h, dtx])

        # generate IRS points
        self.pIRS = self.IRS()#(M, N, dx)

        ## generate sat trajectory
        self.pR0 = np.array([0, satH-h, drx])
        self.rO = er + satH # height of orbit

        if T is None:
            self.T = np.linspace(*self.elev(10), Tpoints)
            if ~np.any(self.T==0):
                # make sure 0 is included to have the maximum present
                self.T = np.append(self.T, 0)
                self.T.sort()
        else:
            self.T = T

        self.pRx = np.asarray([self.pR0 + self.sat_xdir(t) for t in self.T])
        self.pRx = self.pRx @ Rx(rot).T

        self.tau0 = np.linalg.norm(self.pRx-self.pt,axis=1) / c0

        self.tau1 = np.empty((self.pIRS.shape[0],self.pRx.shape[0]))
        for i in range(self.pIRS.shape[0]):
            self.tau1[i] = (np.linalg.norm(self.pRx-self.pIRS[i],axis=-1) + np.linalg.norm(self.pIRS[i]-self.pt))/c0

        self.A0 = self._A0()
        self.Amn = self._Amn()


    # Antenna gains
    def theta(self, p1, p2): # polar
        if p1.ndim == 2 or p2.ndim == 2:
            p1 = _reshape(p1)
            p2 = _reshape(p2)
        elif not ( p1.ndim == 1 and p2.ndim == 1):
            raise ValueError("only 1D and 2D arrays supported")

        return np.arccos((p2[2]-p1[2]) / np.linalg.norm(p2-p1, axis=0))


    def phi(self, p1, p2): # azimuth
        if p1.ndim == 2 or p2.ndim == 2:
            p1 = _reshape(p1)
            p2 = _reshape(p2)
        elif not ( p1.ndim == 1 and p2.ndim == 1):
            raise ValueError("only 1D and 2D arrays supported")

        return np.arctan2(p2[1] - p1[1], p2[0] - p1[0])


    def G(self, p): # IRS
        th = self.theta(self.pIRS, p)

        if (self.isoIRS):
            ret = np.ones_like(th)
            ret[th > np.pi/2] = 0
            return ret
        else:
            ph = self.phi(self.pIRS, p)

            A = 4 * np.pi * self.dx * self.dy * (self.fc / c0)**2
            
            ang = np.cos(th)
            ang[ang<0] = 0

            return A * ang


    def GT(self, p):
        th = self.theta(self.pt, p)
        ph = self.phi(self.pt, p)
        return np.ones_like(th)

        #if p.ndim == 2:
        #    p1 = _reshape(pt)
        #    p = _reshape(p)
        #elif p.ndim > 2:
        #    raise ValueError("only 1D and 2D arrays supported")
        #else:
        #    p1 = pt
        #return np.pi * (p[1]-p1[1]) / np.linalg.norm(p-p1, axis=0)


    def GR(self, pR, pD):
        th = self.theta(pR, pD)
        #ph = phi(pR, pD)
        return np.ones_like(th)


    # generate IRS element centers
    def IRS(self):
        def interval(m):
            mi = (self.M+1)%2 -np.floor(m/2)
            ma = np.floor(m/2)

            assert(mi.is_integer() and ma.is_integer())
            return (int(mi), int(ma))

        Mmin, Mmax = interval(self.M)
        Nmin, Nmax = interval(self.N)

        return np.asarray([(m*self.dx - .5 * self.dx * ((self.M+1)%2), n*self.dy - .5 * self.dy * ((self.N+1)%2), 0) for m, n in it.product(range(Mmin, Mmax+1), range(Nmin, Nmax+1))])


    # satellite trajectory (moves in positive x direction)
    def sat_xdir(self, t):
        tmp = np.sqrt(GM/self.rO**3)
        return [self.rO*np.sin(tmp*t), self.rO*(np.cos(tmp*t)-1), 0]


    def sat_xdir_diff(self, t):
        v = np.sqrt(GM/self.rO)
        return [v*np.cos(v/self.rO * t), -v*(np.sin(v/self.rO * t)-1), 0]


    # amplitude gains
    def _A0(self):
        a0 = c0/(4*np.pi*self.fc)

        GTR = self.GT(self.pRx)
        GRT = self.GR(self.pRx, self.pt)

        return a0 * np.sqrt(GTR * GRT) / np.linalg.norm(self.pRx-self.pt, axis=1)


    def _Amn(self):
        GTmn = self.GT(self.pIRS)
        GmnT = self.G(self.pt)

        airs = np.sqrt(self.mu) * (c0/(4*np.pi*self.fc))**2 * np.sqrt(GTmn*GmnT) / np.linalg.norm(self.pIRS-self.pt, axis=1)

        Amn = np.empty((self.pRx.shape[0],))
        for t in range(len(Amn)):
            GmnR = self.G(self.pRx[t])
            GRmn = self.GR(self.pRx[t], self.pIRS)

            Amn[t]  = np.sum(airs * np.sqrt(GmnR * GRmn) / np.linalg.norm(self.pIRS - self.pRx[t], axis=1))

        return Amn


    def zeroh(self):
        tmp = np.linalg.solve(Rx(self.rotation), self.pt)
        return np.arccos((tmp[1]-self.pR0[1]+self.rO)/self.rO) * self.rO/np.sqrt(GM/self.rO)

    def elevation(self, t):
        tmpp = self.pR0 - np.linalg.solve(Rx(self.rotation), self.pt)
        p = tmpp + self.sat_xdir(t)
        return np.pi/2 - np.arccos(p[1]/np.linalg.norm(p))

    def elev(self, a):
        minel = np.pi/2 - a*np.pi/180
        tmpp = self.pR0 - np.linalg.solve(Rx(self.rotation), self.pt)

        def elevang(t):
            p = tmpp + self.sat_xdir(t)
            return minel - np.arccos(p[1]/np.linalg.norm(p))

        zh = self.zeroh()

        return (scipy.optimize.brentq(elevang, -zh, 0), scipy.optimize.brentq(elevang, zh, 0))

    def elevations(self):
        tmp = self.pRx - self.pt
        return np.pi/2 - np.arccos(tmp @ Rx(self.rotation) @ np.array([0,1,0]) / np.linalg.norm(tmp,axis=-1))

    def phaseshifts(self):
        tmp = fc*(self.tau0 - self.tau1)

        return 2*np.pi * (tmp + np.ceil(-tmp))

    def phaseshifts2(self):
        return 2*np.pi * np.mod(fc*(self.tau0 - self.tau1), 1)

    def delaySpread(self):
        Td = np.ceil(self.fc*(self.tau1-self.tau0))/self.fc
        return np.max(Td,axis=0)


class rndSatIRS(SatIRS):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _A0(self):
        A0 = super()._A0()

        return A0 * np.exp(-1j*2*np.pi*self.fc*self.tau0)

    def _Amn(self):
        GTmn = self.GT(self.pIRS)
        GmnT = self.G(self.pt)

        airs = np.sqrt(self.mu) * (c0/(4*np.pi*self.fc))**2 * np.sqrt(GTmn*GmnT) / np.linalg.norm(self.pIRS-self.pt, axis=1)

        phi = np.random.uniform(high=2*np.pi, size=(self.tau1.shape[0],1))
        phases = np.exp(-1j*2*np.pi*self.fc*self.tau1-1j*phi).T
        Amn = np.empty((self.pRx.shape[0],), dtype=np.complex)
        for t in range(len(Amn)):
            GmnR = self.G(self.pRx[t])
            GRmn = self.GR(self.pRx[t], self.pIRS)

            Amn[t]  = np.sum(airs * np.sqrt(GmnR * GRmn) / np.linalg.norm(self.pIRS - self.pRx[t], axis=1) * phases[t])

        return Amn


class Specular(SatIRS):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _A0(self):
        A0 = super()._A0()

        return A0 * np.exp(-1j*2*np.pi*self.fc*self.tau0)

    def _Amn(self):
        GTmn = self.GT(self.pIRS)
        GmnT = self.G(self.pt)

        airs = np.sqrt(self.mu) * (c0/(4*np.pi*self.fc))**2 * np.sqrt(GTmn*GmnT) / np.linalg.norm(self.pIRS-self.pt, axis=1)

        k = 2 * np.pi * self.fc / (c0 * np.linalg.norm(self.pIRS-self.pt, axis=1))
        k = k[:,np.newaxis] * ((self.pIRS-self.pt) * [-1,-1,1])

        Amn = np.empty((self.pRx.shape[0],), dtype=np.complex)
        for t in range(len(Amn)):
            phases = np.exp(-1j*(np.sum((self.pIRS - self.pRx[t]) * k, axis=1) + 2*np.pi*self.fc*self.tau0[t]))

            GmnR = self.G(self.pRx[t])
            GRmn = self.GR(self.pRx[t], self.pIRS)

            Amn[t]  = np.sum(airs * np.sqrt(GmnR * GRmn) / np.linalg.norm(self.pIRS - self.pRx[t], axis=1) * phases)

        return Amn

if __name__=="__main__":
    import time
    tic = time.time()

    MC = 10 # average over MC random inits
    #MC = 2
    Tpoints = 500
    #Tpoints = 10
    fc = 2e9
    lc = c0 / fc
    dx = dy = lc / 5

    drx = 1000
    dtx = 1000

    # US Billboard (largest Bulletin: https://dashtwo.com/blog/what-are-billboard-dimensions/)
    #N = int(np.round(14 * .3048 / dy))
    #M = int(np.round(48 * .3048 / dx))
    #spectacular
    N = int(np.round(20 * .3048 / dy))
    M = int(np.round(60 * .3048 / dx))

    N = int(np.round(40 * .3048 / dy))
    # general
    #M = N = 1000

    iso = SatIRS(M, N, dx, dy, rot = 0, fc = fc, isoIRS = True, drx = drx, dtx = dtx, Tpoints=Tpoints)
    planar = SatIRS(M, N, dx, dy, rot = 0, fc = fc, drx = drx, dtx = dtx, Tpoints=Tpoints)
    tilt = SatIRS(M, N, dx, dy, rot = np.pi/4, fc = fc, drx = drx, dtx = dtx, Tpoints=Tpoints)

    specular = Specular(M, N, dx, dy, rot = 0, fc = fc, drx = drx, dtx = dtx, Tpoints=Tpoints)
    speculartilt = Specular(M, N, dx, dy, rot = np.pi/4, fc = fc, drx = drx, dtx = dtx, Tpoints=Tpoints)

    planarRndGain = np.empty((MC,len(planar.T)))
    tiltRndGain = np.empty((MC,len(tilt.T)))
    for i in range(MC):
        planarRnd = rndSatIRS(M, N, dx, dy, rot = 0, fc = fc, drx = drx, dtx = dtx, Tpoints=Tpoints)
        tiltRnd = rndSatIRS(M, N, dx, dy, rot = np.pi/4, fc = fc, drx = drx, dtx = dtx, Tpoints=Tpoints)

        planarRndGain[i] = np.abs(planarRnd.A0+planarRnd.Amn)
        tiltRndGain[i] = np.abs(tiltRnd.A0+tiltRnd.Amn)
    planarRndGainAvg = np.mean(planarRndGain, axis=0)
    tiltRndGainAvg = np.mean(tiltRndGain, axis=0)

    #pd.DataFrame(data = 20*np.log10(iso.A0), index = iso.T, columns = 'isotropic IRS')

    #g = ((tmp.A0 + tmp.Amn)**2 / tmp.A0**2 - 1) * 100
    #print("min: {:.4g}%\nmax: {:.4g}%".format(np.min(g), np.max(g)))

    # sanity check
    assert(np.allclose(iso.T, planar.T))
    assert(np.allclose(iso.T, tilt.T))
    assert(np.allclose(iso.elevations(), planar.elevations()))
    assert(np.allclose(iso.elevations(), tilt.elevations()))

    # create data
    base = pd.DataFrame({'Time': iso.T, 'Elevation': iso.elevations()*180/np.pi})
    #base = base.set_index('Time')
    data = pd.DataFrame({
        'no IRS': iso.A0,
        'isotropic IRS': iso.A0+iso.Amn,
        'planar IRS': planar.A0+planar.Amn,
        'tilt IRS': tilt.A0+tilt.Amn,
        'tilt IRS no LOS': tilt.Amn,
        'random': planarRndGainAvg,
        'random tilt': tiltRndGainAvg,
        'specular': np.abs(specular.A0 + specular.Amn),
        'specular tilt': np.abs(speculartilt.A0 + speculartilt.Amn),
        })
    data = 20*np.log10(data)

    channel = base.join(data)
    gain = base.join(data[['isotropic IRS','planar IRS','tilt IRS', 'tilt IRS no LOS', 'random', 'random tilt', 'specular', 'specular tilt']].apply(lambda x: x-data['no IRS']))

    channel = channel.set_index('Time')
    gain = gain.set_index('Time')
    
    # print some values
    maxmin = pd.DataFrame({'min gain': gain.min(), 'max gain': gain.max()})
    tmp = (10**(maxmin[2:]/20)-1)*100
    maxmin = maxmin.add_suffix(' dB').join(tmp.add_suffix(' %'))
    print(maxmin)

    # generate time to elevation table
    ang = [10, 25, 50, 90]
    elev = pd.DataFrame()
    for a in ang:
        t = np.around(iso.elev(a), decimals = 3)

        if t[0] == t[1]:
            t = np.delete(t, 1)
        else:
            assert(t[0] == -t[1])

        elev = elev.append(pd.DataFrame([a]*len(t),index=t))
    elev = elev.sort_index()
    elev = elev.rename(columns={0: 'Elevation'})

    # store data
    channel.to_csv('../channel.dat')
    gain.to_csv('../gain.dat')
    elev.to_csv('../elevation.dat')

    # plot data
    fig = plt.figure()

    ax = fig.add_subplot(121)
    ax.set_title('Channel Gain')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Gain [dB]')
    channel.plot(ax=ax, y=data.columns, color=['C' + str(d) for d in range(data.shape[1])])

    ax = fig.add_subplot(122)
    ax.set_title('IRS Gain over Baseline')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Gain [dB]')
    gain.plot(ax=ax, y=data.columns.drop('no IRS'), color=['C' + str(d) for d in range(data.shape[1])][1:], legend = None)

    toc = time.time()
    plt.show()


#tmp = np.linalg.norm(pRx-pt, axis=1)
#tau0 = tmp / c0
#Ds = np.sum(((pRx - pt).T/tmp).T * np.asarray([sat_xdir_diff(t) for t in T]), axis=1) * fc # not so sure about this...

#Airs = np.asarray([np.sum(airs / np.linalg.norm(pIRS-pRx[1], axis=1)) for t in range(len(pRx))])

## DISPLAY

# plot IRS as surface
#def plotIRS(ax, points, *args, **kwargs):
#    pIRS = points.T
#    X, Y = np.meshgrid([np.min(pIRS[0]), np.max(pIRS[0])], [np.min(pIRS[0]), np.max(pIRS[0])])
#    Z = np.zeros_like(X)
#
#    ax.plot_surface(X, Z, Y, *args, **kwargs)
#    
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection = '3d')
#
## exchange y and z axis!
#ax.set_xlabel('x')
#ax.set_ylabel('z')
#ax.set_zlabel('y')
#
#
#ax.scatter(pt[0], pt[2], pt[1], label='Tx')
#ax.scatter(0, 0, 0, label='IRS', color='green')
#plotIRS(ax, pIRS, color='green')
#
#ax.scatter(pR0[0], pR0[2], pR0[1], label='Rx')
#ax.plot(pRx[:,0], pRx[:,2], pRx[:,1])

#ax.legend()

#plt.show()
