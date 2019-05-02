import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


#user inputs

n = int(input("Please enter lattice dimension, n x n: "))
dynamic = str(input("Please choose GLAUBER or KAWASAKI dynamics: "))
dynamic = dynamic.upper()

T_Lower = float(input("Please enter a lower temperature value: "))
T_Upper = float(input("Please enter an upper temperature value: "))
T_animate = float(input("Please enter a temperature to animate: "))

Animate = str(input("Type YES for animaiton, NO for no animaiton: "))
Animate = Animate.upper()

passes = int(input("Please enter the number of sweeps: "))

n_mes = int(input("Please enter the number of measurements per temperature: "))

#==============================================================================

N = n*n
E_list = []
M_list = []
Cv_1 = []
Cv_sigma = []
X_1 = []
X_sigma = []
M_sigma = []
E_sigma = []
nt = 30
mc_step = 1000
T =  np.linspace(T_Lower,T_Upper,nt)
 

#lattice       
def initialState(n):
    down = -1
    
    i_matrix = np.random.randint(0,1+1,size=(n,n))
    
    for i in range(n):
        for j in range(n):
            if i_matrix[i][j] == 0:
                i_matrix[i][j] = down
    return i_matrix

   
    
#glauber dynamics function
def glauberDynamic(config,te,mc_step,n):

        if dynamic == "GLAUBER":

            for i in range(mc_step):
                deltaE = 0
                a = int(np.random.uniform(0, n))
                b = int(np.random.uniform(0, n))
                s =  config[a, b]
                nb = config[(a+1)%(n-1),b] + config[a,(b+1)%(n-1)] + config[(a-1)%(n-1),b]
                + config[a,(b-1)%(n-1)]
                deltaE = 2*s*nb
                
                if deltaE <= 0:
                    s *= -1
                elif np.random.rand() < np.exp(-(1*deltaE)/te):
                    s *= -1
                config[a, b] = s
        return config


#kawasaki dynamics function
def kawasakiDynamic(config,te,mc_step,n):

        for i in range(mc_step):

            a1 = int(np.random.randint(0, n))
            b1 = int(np.random.randint(0, n))
            a2 = int(np.random.randint(0, n))
            b2 = int(np.random.randint(0, n))

            s1 =  config[a1,b1]
            s2 =  config[a2,b2]
            
            if s1 != s2:
                nb1 = config[(a1+1)%(n-1),b1] + config[a1,(b1+1)%(n-1)] 
                + config[(a1-1)%(n-1),b1] + config[a1,(b1-1)%(n-1)]
                nb2 = config[(a2+1)%(n-1),b2] + config[a2,(b2+1)%(n-1)] 
                + config[(a2-1)%(n-1),b2] + config[a2,(b2-1)%(n-1)] 
                deltaE = 2*s1*nb1* + 2*s2*nb2

                
                if a1+1==a2 and b1==b2: deltaE+=4
                elif a1-1==a2 and b1==b2: deltaE+=4
                elif a1==a2 and b1+1==b2: deltaE+=4
                elif a1==a2 and b1-1==b2: deltaE+=4
                
                
                if deltaE <= 0:
                    config[a1, b1] = s2
                    config[a2, b2] = s1
                elif np.random.rand() < np.exp(-(1*deltaE)/te):
                    config[a1, b1] = s2
                    config[a2, b2] = s1

        return config
    
#lattice energy calculation 
def latticeEnergy(config,n):
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%(n-1),j] + config[i,(j+1)%(n-1)] 
            + config[(i-1)%(n-1),j] + config[i,(j-1)%(n-1)]
            energy += -1*S*nb
    e2 = energy**2

    return energy/4,e2/16

#susceptibility calculation
def susCalc(a1,a2,te,N):    
    avg_m = np.mean(a1)
    avg_m2 = np.mean(a2)
    X = (avg_m2 - avg_m*avg_m)/(te*N)
    return avg_m, X

#magnetization calculation
def magCalc(config,sq):
    if sq == False:
        m = np.sum(config)
    elif sq == True:
        m = (np.sum(config))**2
    return np.abs(m)

#error susceptibility using bootstrap
def errorSus(sx,num,samp):
    ar = np.zeros(num-1)
    sam = np.zeros(samp)
    for k in range(samp):
        for j in range(len(ar)):
            x = np.random.randint(0,high=num)
            ar[j] = sx[x]
        avg = np.mean(ar)
        sam[k] = avg
    cmean = np.mean(sam)**2
    c2mean = np.mean(sam*sam)
    erc = np.sqrt((c2mean-cmean))
    return erc

#error magnetization
def errorMag(m,num):
    er1 = np.std(m)/(len(m))**0.5
    return er1

#animation function
def Simulation(*args):
    if dynamic == "GLAUBER":
        im.set_data(glauberDynamic(ani_config,T_animate,100,n))
    elif dynamic == "KAWASAKI":
        im.set_data(kawasakiDynamic(ani_config,T_animate,1000,n))        
    return im,




#==============================================================================

#main function
if Animate != "YES":
    z = 0
    config = initialState(n)
    for i in T:
        temp = i
        sweeps = 0
        m_list = []
        m2_list = []
        e_list = []
        e2_list = []
        Cv_list = np.zeros(n_mes)
        e_Avg = np.zeros(n_mes)
        m_avg = np.zeros(n_mes)
        X_avg = np.zeros(n_mes)

        for k in range(n_mes):
            for j in range(passes):
                if dynamic == "GLAUBER":

                    glauberDynamic(config,temp,mc_step,n)
                elif dynamic == "KAWASAKI":
                    kawasakiDynamic(config,temp,mc_step,n)


                if sweeps >= 500:
                    if sweeps%10 == 0:
                        m_list.append(magCalc(config,sq=False))
                        m2_list.append(magCalc(config,sq=True))
                        e,e2 = latticeEnergy(config,n)
                        e_list.append(e)
                        e2_list.append(e2)
                sweeps += 1
    
            avg_E, CvN = susCalc(e_list,e2_list,temp,N)
            avg_mag, X = susCalc(m_list,m2_list,temp,N)
            
            e_Avg[k] = avg_E/N
            Cv_list[k] = (CvN/temp)
            X_avg[k] = (X)
            m_avg[k] = avg_mag/N
    
        Cv_sigma.append(errorSus(Cv_list,n_mes,1000))
        X_sigma.append(errorSus(X_avg,n_mes,1000))
        M_sigma.append(errorMag((m_avg),n_mes))
        E_sigma.append(errorMag((e_Avg),n_mes))

    

    
        M_list.append(np.mean(m_list)/N)
        X_1.append(X)
        Cv_1.append(CvN/temp)
        E_list.append(np.mean(e_list)/N)
        z += 1
  
    
#==============================================================================
#animation 
if Animate == "YES":
    ani_config = initialState(n)

    fig = plt.figure()
    im = plt.imshow(kawasakiDynamic(ani_config,T_animate,mc_step,n), animated=True)
    plt.title(T_animate)
    
    for i in range(1000):
        Simulation()
    
    ani = animation.FuncAnimation(fig, Simulation ,interval=25, blit=True
                                  ,frames = 500,repeat = True)
    plt.show()
#==============================================================================



#write data to the output file   
if Animate != "YES":
    
    if dynamic == "GLAUBER":
        d = 'GLAUBER'
    else:
        d = 'KAWASAKI'
    
    output = open(str(d) + '_data.txt','w')
    
    output.write('Energy\t' + 'Sigma_E\t' +'Magnitisation\t'
                 +'Sigma_M\t'+'Heat Capacity\t'+'Sigma_C\t'+'Susceptibilty\t'
                 +'Sigma_X\t\n')
    for i in range(len(T)):
        output.write(str(E_list[i]) + '\t' + str(E_sigma[i]) + '\t' + str(M_list[i]) + '\t' + str(M_sigma[i])
        + '\t' + str(Cv_1[i]) + '\t' + str(Cv_sigma[i]) + '\t' + str(X_1[i]) + '\t' + str(X_sigma[i])
        + '\n' )
        
    output.close()
#==============================================================================

    
#plotting functions    
    
    f = plt.figure(figsize=(15, 10))  
    
    ax =  f.add_subplot(2, 2, 1 );
    plt.title('Energy vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('<E>')
    ax.errorbar(T,E_list,yerr = E_sigma,fmt='-', ecolor='orangered',color='steelblue',
                capsize=2)   
    ax =  f.add_subplot(2, 2, 2 );
    plt.title('Specific Heat Capacity vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('<cV/N>')
    ax.errorbar(T,Cv_1,yerr=Cv_sigma,fmt='-', ecolor='orangered',color='steelblue',
                 capsize=2)
    ax =  f.add_subplot(2, 2, 3 );
    plt.title('Average magnitisation vs temperature')
    plt.xlabel('Temperature')
    plt.ylabel('<M>')
    ax.errorbar(T,M_list,yerr = M_sigma,fmt='-', ecolor='orangered',color='steelblue',
                capsize=2)
    
    ax =  f.add_subplot(2, 2, 4 );
    plt.title('Magnetic Susceptibilty vs Temperature')
    plt.xlabel('Temperature')
    plt.ylabel('<X>')
    ax.errorbar(T,X_1,yerr=X_sigma,fmt='-', ecolor='orangered',color='steelblue',
                capsize=2)
    
    #plt.savefig(str(d) +'_all.png')
    plt.show()

