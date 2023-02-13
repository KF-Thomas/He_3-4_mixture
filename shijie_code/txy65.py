import csv
import matplotlib.pyplot as plt
import os
import multiprocessing as mp
import datetime
import numpy as np
import pyximport; pyximport.install() #can direcly use '.so' file instead.
#import ctry2

#----------------
#Initialisation
#Note that if the bin size here changed, the bin size in ctry2 should also be changed.
#[delta_t,delta_x,delta_y]=[0.002,0.1,0.1]
[delta_t,delta_x,delta_y]=[0.002,0.00013,0.00056] #bin size.
#[delta_t,delta_x,delta_y]=[0.00005,0.0005,0.0005]
#path='/Users/shijieli/Desktop/project1/csvfiles2/' #folder path
#path='/Users/shijieli/Desktop/project1/csvfile_fermion_2/'
#path='/Users/shijieli/Desktop/project1/tf/20220325_TF_dep_trial_1/csvfile/20220325_he_3_4_mix_evap_840_detune_765/'
path='Z:/EXPERIMENT-DATA/2022_degenerate_He3_and_He4_mixture/he3_he4_mixture/20220325_TF_dep_trial_1/20220325_he_3_4_mix_evap_840_detune_765/'
bins=15
hist_upper=delta_t
[t_min,t_max]=[0.4,0.44]
#----------------
def takeSize(i):
    if path[-1] == '/':
        return os.path.getsize(path+i)
    else:
        return os.path.getsize(path+'/'+i)

def task0(file):
    with open(path+file) as csvfile:
        reader = csv.reader(csvfile)
        table = [ row for row in reader]
        txy_temp=np.loadtxt(path+file,delimiter=',')#For csv files. For other types of files just change the delimiter.
        txy=txy_temp[np.where((txy_temp[:,0]>t_min) & (txy_temp[:,0]<t_max))]
    pairs=ctry2.task(txy)[0]
    return [txy,pairs]

def plot_G(a):
    a[0][0]=a[0][0]*(1+1e-07/(delta_t/bins))
    plt.cla()
    plt.ylim(0.9*min(a[0]),1.1*max(a[0]))
    #plt.hist(pairs,bins)
    plt.scatter(a[1][1:],a[0])
    #plt.scatter(a[1][1:],a[0])
    plt.xlabel(r'$\delta \tau(s)$')
    plt.ylabel(r'$G^{(2)}(\tau)$')
    plt.grid()
    plt.show()

def plot_g(a,b):
    g=[a[0][i]*N/b[i] for i in range(len(a[0]))]
    t=[i*1000 for i in a[1][:-1]]
    plt.cla()
    plt.xlabel(r'$\delta T(ms)$')
    plt.ylabel('g')
    #plt.errorbar(t,g,yerr=error_0,fmt='o',ms=4,capsize=3)
    plt.scatter(t,a)
    plt.grid()
    plt.show()
    
def pairs_index(x):
    x=int(x)
    a=[]
    for i in range(x-1):
        for m in range(i+1,x):
            a.extend([[i,m]])
    return np.array(a)

def nomalization_profile(txy_pool):
    pairs_length=[pow(file_length[i[0]]*file_length[i[1]],1.85/2) for i in indexs]
    N=sum(pairs_length)/sum([pow(i,1.85) for i in file_length])
    return N

def write_file(list0,name):
    with open(path+'data/'+name+'.txt','w') as csvfile:
        writer=csv.writer(csvfile)
        for i in list0:
            aa=writer.writerow([str(i)])
    return

def plot_errorbar(x,y,error):
    plt.errorbar(x,y,yerr=error,fmt='o',ms=4,capsize=3)
    plt.xlabel(r'$\delta t(s)$')
    plt.ylabel(r'$g(\delta t)$')
    plt.grid()
    plt.show()

def read_file(file):
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        table = [ float(row[0]) for row in reader]
    return table

def task_2_1(hists):
    aa_power=[aa[n]/file_length_power[n] for n in range(len(aa))]
    error=[np.std([i[n] for i in aa_power]) for n in range(len(hists[0]))]
    nor=[sum(hist_G[0][:i+1]) for i in range(bins)]
    error_G=[error[i]*np.sqrt(len(files2))/nor[i] for i in range(len(nor))]
    nor1=[sum(hist_g[:i+1]) for i in range(bins)]
    y=[nor[i]/nor1[i]*N for i in range(bins)]
    error_0=[y[i]*error_G[i] for i in range(bins)]
    return [error_0,y]

def task_2_2():
    aa_power=[aa[n]/file_length_power[n] for n in range(len(aa))]
    error=[np.std([i[n] for i in aa_power]) for n in range(bins)]
    error_G=[error[i]*np.sqrt(len(files2))/hist_G[0][i] for i in range(bins)]
    error_0=[y[i]*error_G[i] for i in range(bins)]
    return error_0

def task_2_2_test(aa,file_length_power):
    aa_power=[aa[n]/file_length_power[n] for n in range(len(aa))]
    aa_av=sum(aa_power)/len(aa)
    error=[np.std([i[n] for i in aa_power]) for n in range(bins)]
    error_G=np.array(error)/aa_av/np.sqrt(len(aa))
    return error_G

def task_2_3():
    x=[1+0.01*i for i in range(310)]
    y=[]
    for n in x:
        sum_0=sum([pow(i,n) for i in file_length])
        file_length_power=[pow(i,n)*len(files2)/sum_0 for i in file_length]
        y.append(task_2_2_test(aa,file_length_power)[0])
    plt.plot(x,y)
    plt.grid()
    plt.ylabel(r'$\delta \overline{G^{(2)}_k}/\overline{G^{(2)}_k}$')
    plt.xlabel('n')
    plt.show()

def tlist(bins0,tmin,tmax):#400
    t_hist=np.zeros(bins0)
    for i in txy_pool:
        t_hist+=np.histogram(i[:,0],bins0,(tmin,tmax))[0]
    time=np.histogram(txy_pool[0][:,0],bins0,(tmin,tmax))[1][:-1]
    plt.plot(time,t_hist/len(files2)*10)
    plt.ylabel(r'Flux(kHz)')
    plt.xlabel(r'Time(s)')
    plt.grid()
    plt.show()
    return list(t_hist/len(files2))

def xlist(bins0):#400
    t_hist=np.zeros(bins0)
    for i in txy_pool:
        t_hist+=np.histogram(i[:,1],bins0,(-0.03,0.03))[0]
    time=np.histogram(txy_pool[1][:,1],bins0,(-0.03,0.03))[1][:-1]
    plt.plot(time,t_hist/len(files2)*10)
    plt.ylabel(r'Flux(kHz)')
    plt.xlabel(r'Time(s)')
    plt.grid()
    plt.show()
    #return list(t_hist/len(files2)*10)
    return list(time)

def hist(bins):
    hist_g=np.histogram(pairs_g,bins,(0,hist_upper))
    hist_G=np.histogram(pairs_G,bins,(0,hist_upper))
    N=nomalization_profile(txy_pool)
    y=[hist_G[0][i]*N/hist_g[0][i] for i in range(bins)]
    #[error_0,y]=task_2_1(aa)
    error_0=task_2_2(hist_g)
    plot_g(hist_G,hist_g)#506

if __name__ == '__main__':
    files2=[i for i in os.listdir(path) if i[:2]=='d_']# and os.path.getsize(path+i)>10000]#read csv files like 'data_1.csv'
    files2.sort(key=takeSize)
    files2.reverse()
    #indexs=pairs_index(len(files2))
    #------------------
    #If there are 1,600 files.
    files2=files2[0:1600]
    indexs0=pairs_index(100)
    indexs=indexs0
    for i in range(1,16):
        indexs=np.concatenate((indexs,indexs0+i*100),axis=0)
    #------------------
    indexs=[[i[1],i[0]] for i in indexs]+[i for i in indexs]
    indexs=[]
    num_cores = int(mp.cpu_count())
    print("Cores:" + str(num_cores)+'\nTask1:' +str(len(files2)))
    pool = mp.Pool(num_cores) #Cheack the number of cores.
    txy_pool=[]
    results=[pool.apply_async(task0, args=[names]) for names in files2]
    pairs_G=[]
    aa=[]
    hist_G=[]
    start_t = datetime.datetime.now()
    for p in results: #Parallel calculating for G2
        a=p.get()
        txy_pool.append(a[0])
        aa.append(np.histogram(a[1],bins,(0,hist_upper))[0])
        pairs_G.extend(a[1])
    file_length=[i.shape[0] for i in txy_pool]
    sum_0=sum([pow(i,1.85) for i in file_length])
    file_length_power=[pow(i,1.85)*len(files2)/sum_0 for i in file_length]#1.85
    print('Task2:'+str(len(indexs)))
    results=[pool.apply_async(ctry2.task_2, args=[txy_pool[i[0]],txy_pool[i[1]]]) for i in indexs]
    hist_g=np.zeros(bins)
    #pairs_g=read_file(path+'data/'+'small_g_2022-08-24 20:09:36.txt') #read files if data saved.
    for p in results: #Parallel calculating for normalisation profile
        a=p.get()
        hist_g+=a
    end_t = datetime.datetime.now()
    #write_file(pairs_G,'big_G_'+str(end_t)[0:19]) #save files
    #write_file(pairs_g,'small_g_'+str(end_t)[0:19])
    elapsed_sec = (end_t - start_t).total_seconds()
    pool.close()
    print("Calculation time: " + "{:.2f}".format(elapsed_sec) + " seconds")
    #hist_g=np.histogram(pairs_g,bins,(0,hist_upper))
    hist_G=np.histogram(pairs_G,bins,(0,hist_upper))
    N=nomalization_profile(txy_pool)#N=nomalization_profile(txy_pool) This command is important
    y=[hist_G[0][i]*N/hist_g[i] for i in range(bins)]#2d
    #[error_0,y]=task_2_1(aa)#For x or y axis
    error_0=task_2_2()#2d
    plot_errorbar(hist_G[1][1:],y,error_0)
    #tlist(1000,0.4,0.5) #TOF (bins,tmin,tmax)
