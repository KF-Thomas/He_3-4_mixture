import csv
import matplotlib.pyplot as plt
import os
import multiprocessing as mp
import datetime
import numpy as np
import pyximport; pyximport.install()
import ctry2

#[delta_t,delta_x,delta_y]=[0.002,0.1,0.1]
#[delta_t,delta_x,delta_y]=[0.004,0.000042,0.00042]
[delta_t,delta_x,delta_y]=[0.006,0.00013,0.00056]
#[delta_t,delta_x,delta_y]=[0.00005,0.0005,0.0005]
#path='/Users/shijieli/Desktop/project1/csvfile_fermion/'
path='/Users/shijieli/Desktop/project1/csvfile_fermion_2/'
bins=15
hist_upper=delta_t
[t_min,t_max]=[0.4,0.44]#3.2-3.3 0.4-0.44

def takeSize(i):
    return os.path.getsize(path+i)

def task0(file):
    with open(path+file) as csvfile:
        reader = csv.reader(csvfile)
        table = [ row for row in reader]
        #txy=[[float(i[0]),float(i[1]),float(i[2])] for i in table]
        txy_temp=np.loadtxt(path+file,delimiter=',')
        txy=txy_temp[np.where((txy_temp[:,0]>t_min) & (txy_temp[:,0]<t_max))]
    [pairs,xy]=ctry2.task_4(txy)#task
    return [txy,pairs,xy]

def task(txy1):
    pair=[]
    length=len(txy1)
    [tmin,tmax]=[0,1]
    xy=np.array([[],[]])
    x=np.array([])
    for n in range(length):
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
                temp=np.array([x,x-dt])
                xy=np.concatenate((xy,temp),axis=1)
                x=np.concatenate((x,np.array([dt])),axis=0)
    return pair

def task_2(txy1,txy2):
    pair=[]
    length_txy2=len(txy2)
    length_txy1=len(txy1)
    [tmin,tmax]=[0,1]
    for n in range(length_txy1):
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
                pair+=[abs(txy2[i][0]-txy1[n][0])]
        
    return pair

def plot_G(a):
    #a[0][0]=a[0][0]*(1+1e-07/(delta_t/bins))
    plt.cla()
    plt.ylim(0.9*min(a[0]),1.1*max(a[0]))
    #plt.hist(pairs,bins)
    plt.scatter(a[1][1:],a[0])
    plt.xlabel(r'$\delta T(s)$')
    plt.ylabel('number of pairs')
    plt.grid()
    plt.show()

def plot_g(a,b):
    #g=[a[0][i]*N/(b[0][i]+a[0][i]) for i in range(len(a[0]))]
    g=[a[0][i]*N/b[i] for i in range(len(a[0]))]
    t=[i*1000 for i in a[1][:-1]]
    plt.cla()
    plt.xlabel(r'$\delta T(ms)$')
    plt.ylabel('g')
    #plt.errorbar(t,g,yerr=error_0,fmt='o',ms=4,capsize=3)
    plt.scatter(t,g)
    plt.grid()
    plt.show()
    
def pairs_index(x):
    x=int(x)
    a=[]
    for i in range(x-1):
        for m in range(i+1,x):
            a.extend([[i,m]])
    return a

def pairs_index_3d(x):
    x=int(x)
    a=[]
    for i in range(x-2):
        for m in range(i+1,x-1):
            for n in range(m+1,x):
                a.append([i,m,n])
    return np.array(a)

def pairs_index_4d(x):
    x=int(x)
    a=[]
    for i in range(x-3):
        for m in range(i+1,x-2):
            for n in range(m+1,x-1):
                for p in range(n+1,x):
                    a.append([i,m,n,p])
    return np.array(a)

def nomalization_profile():
    pairs_length=[file_length[i[0]]*file_length[i[1]] for i in indexs]
    N=sum(pairs_length)/sum([i**2 for i in file_length])
    return N

def nomalization_profile_3d():
    pairs_length=[file_length[i[0]]*file_length[i[1]]*file_length[i[2]] for i in indexs]
    N=sum(pairs_length)/sum([i**3 for i in file_length])
    return N

def nomalization_profile_4d():
    pairs_length=[pow(file_length[i[0]]*file_length[i[1]]*file_length[i[2]]*file_length[i[3]],1/4) for i in indexs]
    N=sum(pairs_length)/sum([i**1 for i in file_length])
    return N

def write_file(list0,name):
    with open(path+'data/'+name+'.txt','w') as csvfile:
        writer=csv.writer(csvfile)
        for i in list0:
            aa=writer.writerow([str(i)])
    return

def plot_errorbar(x,y,error):
    plt.errorbar(x,y,yerr=error,fmt='o',ms=4,capsize=3)
    plt.xlabel(r'$\delta x(m)$')
    plt.ylabel(r'$g(\tau=0)$')
    plt.grid()
    plt.show()

def read_file(file):
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        table = [ float(row[0]) for row in reader]
    return table

def task_2_1(hists):
    aa_power=[aa[n]/file_length_power[n] for n in range(len(aa))]
    #aa_power=aa
    error=[np.std([i[n] for i in aa_power]) for n in range(len(hists[0]))]
    nor=[sum(hist_G[0][:i+1]) for i in range(bins)]
    error_G=[error[i]*np.sqrt(len(files2))/nor[i] for i in range(len(nor))]
    nor1=[sum(hist_g[0][:i+1]) for i in range(bins)]
    y=[nor[i]/nor1[i]*N for i in range(bins)]
    error_0=[y[i]*error_G[i] for i in range(bins)]
    return [error_0,y]

def task_2_2(hists):
    aa_power=[aa[n]/file_length_power[n] for n in range(len(aa))]
    error=[np.std([i[n] for i in aa_power]) for n in range(len(hists))]
    error_G=[error[i]*np.sqrt(len(files2))/hist_G[0][i] for i in range(bins)]
    error_0=[y[i]*error_G[i] for i in range(bins)]
    return error_0

def plot_G_3D():
    hist_2d,xedges,yedges=np.histogram2d(xy_pool[0],xy_pool[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])
    plt.clf()
    X,Y=np.meshgrid(xedges[:-1],yedges[:-1])
    ax3=plt.axes(projection='3d')
    ax3.plot_surface(X,Y,hist_2d,cmap='rainbow')
    ax3.set_xlabel(r'$\tau 1$')
    ax3.set_ylabel(r'$\tau 2$')
    ax3.set_zlabel(r'$G(\tau1,\tau2)$')
    plt.show()

def plot_g_3D():
    hist_2d,xedges,yedges=np.histogram2d(xy_pool[0],xy_pool[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])
    plt.clf()
    X,Y=np.meshgrid(xedges[:-1],yedges[:-1])
    #Z=hist_2d/xy_pool2*N
    #Z=xy_pool2
    ax3=plt.axes(projection='3d')
    ax3.plot_surface(X,Y,Z,cmap='rainbow')
    ax3.set_xlabel(r'$\tau 1$')
    ax3.set_ylabel(r'$\tau 2$')
    ax3.set_zlabel(r'$g(\tau1,\tau2)$')
    plt.show()

def plot_g_3D_scatter():
    #Z=hist_2d/xy_pool2*N
    #Z=xy_pool2
    hist_2d,xedges,yedges=np.histogram2d(xy_pool[0],xy_pool[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])
    X,Y=np.meshgrid(xedges[:-1],yedges[:-1])
    ax3=plt.axes(projection='3d')
    ax3.scatter(X,Y,Z)
    ax3.set_xlabel(r'$\tau 1$')
    ax3.set_ylabel(r'$\tau 2$')
    ax3.set_zlabel(r'$g(\tau1,\tau2)$')
    plt.show()

def plot_g_3D_1():
    #hist_2d=np.histogram2d(xy_pool[0],xy_pool[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])
    hist_2d,xedges,yedges=np.histogram2d(xy_pool[0],xy_pool[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])
    plt.clf()
    X,Y=np.meshgrid(xedges[:-1],yedges[:-1])
    #hist_g_2d=np.histogram2d(xy_pool2[0],xy_pool2[1],bins,range=[[0,delta_t/2],[0,delta_t/2]])
    #X,Y=np.meshgrid(hist_2d[1][:-1],hist_2d[2][:-1])
    #X=np.array([X[-n-1] for n in range(len(X))])
    Z=hist_2d/xy_pool2*N
    z00=[Z[i][i] for i in range(Z.shape[0])]
    plt.scatter(X[0],z00)
    plt.xlabel(r'$\delta T(s)$')
    plt.ylabel('g(x=y)')
    plt.grid()
    plt.show()

def plot_g_4D_1():
    plt.clf()
    Z=hist_3d/xy_pool2*N
    X=np.histogramdd((xy_pool[0],xy_pool[1],xy_pool[2]),bins,range=[[0,delta_t/3],[0,delta_t/3],[0,delta_t/3]])[1][1][:-1]
    z00=[Z[i][i][i] for i in range(Z.shape[0])]
    plt.scatter(X,z00)
    plt.xlabel(r'$\delta T(s)$')
    plt.ylabel('g(x=y=z)')
    plt.grid()
    plt.show()

def plot_g_4D_00():
    plt.clf()
    z01_power=[z01[n]/file_length_power[n] for n in range(len(z01))]
    z02=sum(z01)
    Z=hist_3d/xy_pool2*N
    z00=[Z[i][i][i] for i in range(Z.shape[0])]
    #z00=[Z[i][i] for i in range(Z.shape[0])]
    error=[np.std([i[n] for i in z01_power]) for n in range(len(z01[0]))]
    error_G=np.array(error)/(sum(z01_power)/len(z01_power))/np.sqrt(len(files2))
    #error_G=[error[i]*np.sqrt(len(files2))/z02[i] for i in range(bins)]
    error_0=[z00[i]*error_G[i] for i in range(bins)]
    #plt.scatter(X[0],z00)
    plot_errorbar(X,z00,error_0)
    plt.xlabel(r'$\delta T(s)$')
    plt.ylabel('g(x=y)')
    plt.grid()
    return error_0

def plot_g_4D_00_test(file_length_power):
    z01_power=[z01[n]/file_length_power[n] for n in range(len(z01))]
    z02=sum(z01)
    Z=hist_3d/xy_pool2*N
    z00=[Z[i][i][i] for i in range(Z.shape[0])]
    #z00=[Z[i][i] for i in range(Z.shape[0])]
    error=[np.std([i[n] for i in z01_power]) for n in range(len(z01[0]))]
    error_G=np.array(error)/(sum(z01_power)/len(z01_power))/np.sqrt(len(files2))
    return error_G[0]

def tlist(bins0,tmin,tmax):
    t_hist=np.zeros(bins0)
    for i in txy_pool:
        t_hist+=np.histogram(i[:,0],bins0,(tmin,tmax))[0]
    time=np.histogram(txy_pool[0][:,0],bins0,(tmin,tmax))[1][:-1]
    plt.plot(time,t_hist)
    plt.ylabel(r'Hits')
    plt.xlabel(r'Time(s)')
    plt.grid()
    plt.show()
    return [time,t_hist]

def hist(bins):
    hist_g=np.histogram(pairs_g,bins,(0,hist_upper))
    hist_G=np.histogram(pairs_G,bins,(0,hist_upper))
    N=nomalization_profile(txy_pool)
    y=[hist_G[0][i]*N/hist_g[0][i] for i in range(bins)]
    #[error_0,y]=task_2_1(aa)
    error_0=task_2_2(hist_g)
    plot_g(hist_G,hist_g)#506

def task_2_3():
    x=[1+0.01*i for i in range(310)]
    y=[]
    for n in x:
        sum_0=sum([pow(i,n) for i in file_length])
        file_length_power=[pow(i,n)*len(files2)/sum_0 for i in file_length]
        y.append(plot_g_4D_00_test(file_length_power))
    plt.plot(x,y)
    plt.show()
    return 

if __name__ == '__main__':
    start_t = datetime.datetime.now()
    #files=os.listdir(path)
    files2=[i for i in os.listdir(path) if i[:5]=='data_' and os.path.getsize(path+i)>100000]#100000
    files2.sort(key=takeSize)
    files2.reverse()
    files2=files2[:1600]
    #indexs=pairs_index_4d(len(files2))
    indexs0=pairs_index_4d(10)
    indexs=indexs0
    for i in range(1,160):
        indexs=np.concatenate((indexs,indexs0+i*10),axis=0)
    #indexs=[]
    indexs=[[i[1],i[0],i[2],i[3]] for i in indexs]+[[i[0],i[2],i[1],i[3]] for i in indexs]+[[i[0],i[1],i[3],i[2]] for i in indexs]+[[i[0],i[1],i[2],i[3]] for i in indexs]
    num_cores = int(mp.cpu_count())
    print("Cores:" + str(num_cores)+'\nTask1:' +str(len(files2)))
    pool = mp.Pool(num_cores)
    txy_pool=[]
    results=[pool.apply_async(task0, args=[names]) for names in files2]
    pairs_G=[]
    aa=[]
    xy_pool=np.array([[],[],[]])
    hist_G=[]
    xy_pool0=[]
    for p in results:
        a=p.get()
        txy_pool.append(a[0])
        aa.append(np.histogram(a[1],bins,(0,hist_upper))[0])
        pairs_G.extend(a[1])
        xy_pool=np.concatenate((xy_pool,a[2]),axis=1)
        xy_pool0.append(np.histogramdd((a[2][0],a[2][1],a[2][2]),bins,range=[[0,delta_t/3],[0,delta_t/3],[0,delta_t/3]])[0])
    file_length=[i.shape[0] for i in txy_pool]
    sum_0=sum([pow(i,2) for i in file_length])
    file_length_power=[pow(i,2)*len(files2)/sum_0 for i in file_length]
    z01=[]
    for n in xy_pool0:
        z01.append(np.array([n[i][i][i] for i in range(bins)]))
    X=np.histogramdd((xy_pool[0],xy_pool[1],xy_pool[2]),bins,range=[[0,delta_t/3],[0,delta_t/3],[0,delta_t/3]])[1][1][:-1]
    print('Task2:'+str(len(indexs)))
    results=[pool.apply_async(ctry2.task_5, args=[txy_pool[i[0]],txy_pool[i[1]],txy_pool[i[2]],txy_pool[i[3]]]) for i in indexs]
    #xy_pool2=np.array([[],[]])
    xy_pool2=np.zeros((bins,bins,bins))
    hist_g=np.zeros(bins)
    #pairs_g=read_file(path+'data/'+'small_g_2022-08-24 20:09:36.txt')
    for p in results:
        a=p.get()
        #hist_g+=a[0]
        xy_pool2+=a
    end_t = datetime.datetime.now()
    #write_file(pairs_G,'big_G_'+str(end_t)[0:19])
    #write_file(pairs_g,'small_g_'+str(end_t)[0:19])
    elapsed_sec = (end_t - start_t).total_seconds()
    pool.close()
    print("Calculation time: " + "{:.2f}".format(elapsed_sec) + " seconds")
    #hist_g=np.histogram(pairs_g,bins,(0,hist_upper))
    hist_G=np.histogram(pairs_G,bins,(0,hist_upper))
    hist_3d=np.histogramdd((xy_pool[0],xy_pool[1],xy_pool[2]),bins,range=[[0,delta_t/3],[0,delta_t/3],[0,delta_t/3]])[0]
    N=nomalization_profile_4d()#N=nomalization_profile(txy_pool) This command is important
    #y=[hist_G[0][i]*N/hist_g[i] for i in range(bins)]#2d
    #[error_0,y]=task_2_1(aa)#
    #error_0=task_2_2(hist_g)#2d
    #plot_g(hist_G,hist_g)#2d
    #plot_G(hist_G)
    #plot_errorbar(hist_G[1][1:],y,error_0)
    #plt.hist2d(xy_pool[0],xy_pool[1],bins=[10,10],range=[[0,0.002],[0,0.002]])
    #plt.show()
