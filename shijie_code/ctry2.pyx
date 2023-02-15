#cython: language_level=3
#cimport numpy as cnp
import numpy as np
cdef double delta_t=0.002,delta_x=0.00013,delta_y=0.00056
#cdef double delta_t=0.002,delta_x=0.00013,delta_y=0.00042
cdef int bins=15
#cdef double delta_t=0.00005,delta_x=0.0005,delta_y=0.0005


def task(double[:,:] txy1): # for G2 or G3
    cdef list pair=list()
    cdef int tmin=0, tmax=1, length
    cdef double dt
    #cdef np.ndarray[np.float64_t, ndim=2] xy
    length=txy1.shape[0]
    xy=np.array([[],[]])
    for n in range(length):
        x=np.array([])
        for i in range(tmin,length):
            if txy1[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length):
            if txy1[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin,tmax):
            if abs(txy1[i][1]-txy1[n][1])<=delta_x and abs(txy1[i][2]-txy1[n][2])<=delta_y:
                dt=abs(txy1[i][0]-txy1[n][0])
                pair.extend([dt])
                temp=np.array([x,dt-x])
                xy=np.concatenate((xy,temp),axis=1)
                x=np.concatenate((x,np.array([dt])),axis=0)
    return [pair,xy] #xy is G3 and pair is G2

def task_2(double[:,:] txy1, double[:,:] txy2): #for normalisation of g2
    cdef list pair=list()
    cdef int length_txy2,length_txy1
    length_txy2=txy2.shape[0]
    length_txy1=txy1.shape[0]
    cdef int tmin=0,tmax=1
    for n in range(length_txy1):
        x=np.array([])
        for i in range(tmin,length_txy2):
            if txy2[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length_txy2):
            if txy2[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin,tmax):
            if abs(txy2[i][1]-txy1[n][1])<=delta_x and abs(txy2[i][2]-txy1[n][2])<=delta_y:
                dt=abs(txy2[i][0]-txy1[n][0])
                pair.extend([dt])
    hist_g=np.histogram(pair,bins,(0,delta_t))[0]
    return hist_g

def task_3(double[:,:] txy1, double[:,:] txy2, double[:,:] txy3): #for normalisation of g3
    cdef int length_txy2,length_txy1,length_txy3
    cdef int tmin=0,tmax=1,tmin2=0,tmax2=1
    length_txy2=txy2.shape[0]
    length_txy1=txy1.shape[0]
    length_txy3=txy3.shape[0]
    xy=np.array([[],[]])
    for n in range(length_txy1):
        x=np.array([])
        for i in range(tmin,length_txy2):
            if txy2[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length_txy2):
            if txy2[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin2,length_txy3):
            if txy3[i][0]>txy1[n][0]+1e-07:
                tmin2=i
                break
        for i in range(tmax2,length_txy3):
            if txy3[i][0]>txy1[n][0]+delta_t:
                tmax2=i
                break
        for i in range(tmin,tmax):
            if abs(txy2[i][1]-txy1[n][1])<=delta_x and abs(txy2[i][2]-txy1[n][2])<=delta_y:
                dt=abs(txy2[i][0]-txy1[n][0])
                x=np.concatenate((x,np.array([dt])),axis=0)
        for i in range(tmin2,tmax2):
            if abs(txy3[i][1]-txy1[n][1])<=delta_x and abs(txy3[i][2]-txy1[n][2])<=delta_y:
                dt=abs(txy3[i][0]-txy1[n][0])
                temp0=x[np.where(x[:]<dt)]
                temp=np.array([temp0,dt-temp0])
                xy=np.concatenate((xy,temp),axis=1)
    hist=np.histogram2d(xy[0],xy[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])[0]
    return hist

def task_4(double[:,:] txy1): # for G4
    cdef list pair=list()
    cdef int tmin=0, tmax=1, length
    cdef double dt
    #cdef np.ndarray[np.float64_t, ndim=2] xy
    length=txy1.shape[0]
    xy=np.array([[],[],[]])
    for n in range(length):
        x=np.array([])
        y=np.array([])
        x0=np.array([])
        for i in range(tmin,length):
            if txy1[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length):
            if txy1[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin,tmax):
            if abs(txy1[i][1]-txy1[n][1])<=delta_x and abs(txy1[i][2]-txy1[n][2])<=delta_y:
                dt=abs(txy1[i][0]-txy1[n][0])
                pair.extend([dt])
                temp=np.array([x0,y-x0,dt-y])
                xy=np.concatenate((xy,temp),axis=1)
                x=np.concatenate((x,np.array([dt])),axis=0)
                x0=np.concatenate((x0,x[:-1]),axis=0)
                y=np.concatenate((y,np.zeros(x[:-1].shape)+dt),axis=0)
    return [pair,xy]

def task_5(double[:,:] txy1, double[:,:] txy2, double[:,:] txy3, double[:,:] txy4): #for normalisation of g4
    cdef int length_txy2,length_txy1,length_txy3, length_txy4
    cdef int tmin=0,tmax=1,tmin2=0,tmax2=1,tmin3=0,tmax3=1
    cdef double dt
    length_txy2=txy2.shape[0]
    length_txy1=txy1.shape[0]
    length_txy3=txy3.shape[0]
    length_txy4=txy4.shape[0]
    xyz=np.array([[],[],[]])
    for n in range(length_txy1):
        x=np.array([])
        xy=np.array([[],[]])
        for i in range(tmin,length_txy2):
            if txy2[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length_txy2):
            if txy2[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin2,length_txy3):
            if txy3[i][0]>txy1[n][0]+1e-07:
                tmin2=i
                break
        for i in range(tmax2,length_txy3):
            if txy3[i][0]>txy1[n][0]+delta_t:
                tmax2=i
                break
        for i in range(tmin3,length_txy4):
            if txy4[i][0]>txy1[n][0]+1e-07:
                tmin3=i
                break
        for i in range(tmax3,length_txy4):
            if txy4[i][0]>txy1[n][0]+delta_t:
                tmax3=i
                break
        for i in range(tmin,tmax):
            if abs(txy2[i][1]-txy1[n][1])<=delta_x and abs(txy2[i][2]-txy1[n][2])<=delta_y:
                dt=txy2[i][0]-txy1[n][0]
                x=np.concatenate((x,np.array([dt])),axis=0)
        for i in range(tmin2,tmax2):
            if abs(txy3[i][1]-txy1[n][1])<=delta_x and abs(txy3[i][2]-txy1[n][2])<=delta_y:
                dt=txy3[i][0]-txy1[n][0]
                temp0=x[np.where(x[:]<dt)]
                temp=np.array([temp0,np.zeros(temp0.shape)+dt])
                xy=np.concatenate((xy,temp),axis=1)
        for i in range(tmin3,tmax3):
            if abs(txy4[i][1]-txy1[n][1])<=delta_x and abs(txy4[i][2]-txy1[n][2])<=delta_y:
                dt=txy4[i][0]-txy1[n][0]
                temp0=xy[:,np.where(xy[1]<dt)[0]]
                temp=np.array([temp0[0],temp0[1]-temp0[0],dt-temp0[1]])
                xyz=np.concatenate((xyz,temp),axis=1)
    hist=np.histogramdd((xyz[0],xyz[1],xyz[2]),bins,range=[[0,delta_t/3],[0,delta_t/3],[0,delta_t/3]])[0]
    return hist

def task_6(double[:,:] txy1): #for G5
    cdef list pair=list()
    cdef int tmin=0, tmax=1, length
    cdef double dt
    #cdef np.ndarray[np.float64_t, ndim=2] xy
    length=txy1.shape[0]
    xyz=np.array([[],[],[],[]])
    for n in range(length):
        x=np.array([])
        y=np.array([])
        x0=np.array([])
        xy=np.array([[],[],[]])
        temp=np.array([[],[],[]])
        for i in range(tmin,length):
            if txy1[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length):
            if txy1[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin,tmax):
            if abs(txy1[i][1]-txy1[n][1])<=delta_x and abs(txy1[i][2]-txy1[n][2])<=delta_y:
                dt=abs(txy1[i][0]-txy1[n][0])
                pair.extend([dt])
                temp1=np.concatenate((xy,np.array([dt-xy[0]-xy[1]-xy[2]])),axis=0)
                xyz=np.concatenate((xyz,temp1),axis=1)
                x=np.concatenate((x,np.array([dt])),axis=0)
                if x.shape[0]>=2:
                    x0=np.concatenate((x0,x[:-2]),axis=0)
                    y=np.concatenate((y,np.zeros(x[:-2].shape)+x[-2]),axis=0)
                    temp=np.array([x0,y-x0,dt-y])
                    xy=np.concatenate((xy,temp),axis=1)
    return [pair,xyz]

def task_7(double[:,:] txy1, double[:,:] txy2, double[:,:] txy3, double[:,:] txy4, double[:,:] txy5): #for normalisation of g5
    cdef int length_txy2,length_txy1,length_txy3, length_txy4, length_txy5
    cdef int tmin=0,tmax=1,tmin2=0,tmax2=1,tmin3=0,tmax3=1,tmin4=0,tmax4=1
    cdef double dt
    length_txy2=txy2.shape[0]
    length_txy1=txy1.shape[0]
    length_txy3=txy3.shape[0]
    length_txy4=txy4.shape[0]
    length_txy5=txy5.shape[0]
    xyzt=np.array([[],[],[],[]])
    for n in range(length_txy1):
        x=np.array([])
        xy=np.array([[],[]])
        xyz=np.array([[],[],[]])
        for i in range(tmin,length_txy2):
            if txy2[i][0]>txy1[n][0]+1e-07:
                tmin=i
                break
        for i in range(tmax,length_txy2):
            if txy2[i][0]>txy1[n][0]+delta_t:
                tmax=i
                break
        for i in range(tmin2,length_txy3):
            if txy3[i][0]>txy1[n][0]+1e-07:
                tmin2=i
                break
        for i in range(tmax2,length_txy3):
            if txy3[i][0]>txy1[n][0]+delta_t:
                tmax2=i
                break
        for i in range(tmin3,length_txy4):
            if txy4[i][0]>txy1[n][0]+1e-07:
                tmin3=i
                break
        for i in range(tmax3,length_txy4):
            if txy4[i][0]>txy1[n][0]+delta_t:
                tmax3=i
                break
        for i in range(tmin4,length_txy5):
            if txy5[i][0]>txy1[n][0]+1e-07:
                tmin4=i
                break
        for i in range(tmax4,length_txy5):
            if txy5[i][0]>txy1[n][0]+delta_t:
                tmax4=i
                break
        for i in range(tmin,tmax):
            if abs(txy2[i][1]-txy1[n][1])<=delta_x and abs(txy2[i][2]-txy1[n][2])<=delta_y:
                dt=txy2[i][0]-txy1[n][0]
                x=np.concatenate((x,np.array([dt])),axis=0)
        for i in range(tmin2,tmax2):
            if abs(txy3[i][1]-txy1[n][1])<=delta_x and abs(txy3[i][2]-txy1[n][2])<=delta_y:
                dt=txy3[i][0]-txy1[n][0]
                temp0=x[np.where(x[:]<dt)]
                temp=np.array([temp0,np.zeros(temp0.shape)+dt])
                xy=np.concatenate((xy,temp),axis=1)
        for i in range(tmin3,tmax3):
            if abs(txy4[i][1]-txy1[n][1])<=delta_x and abs(txy4[i][2]-txy1[n][2])<=delta_y:
                dt=txy4[i][0]-txy1[n][0]
                temp0=xy[:,np.where(xy[1]<dt)[0]]
                temp=np.concatenate((temp0,np.zeros((1,temp0.shape[1]))+dt),axis=0)
                xyz=np.concatenate((xyz,temp),axis=1)
        for i in range(tmin4,tmax4):
            if abs(txy5[i][1]-txy1[n][1])<=delta_x and abs(txy5[i][2]-txy1[n][2])<=delta_y:
                dt=txy5[i][0]-txy1[n][0]
                temp0=xyz[:,np.where(xyz[2]<dt)[0]]
                temp=np.array([temp0[0],temp0[1]-temp0[0],temp0[2]-temp0[1],dt-temp0[2]])
                xyzt=np.concatenate((xyzt,temp),axis=1)
    hist=np.histogramdd((xyzt[0],xyzt[1],xyzt[2],xyzt[3]),bins,range=[[0,delta_t/4],[0,delta_t/4],[0,delta_t/4],[0,delta_t/4]])[0]
    return hist