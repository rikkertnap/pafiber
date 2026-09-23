#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Imports for the whole notebook
import pandas as pd
import glob as gb
import os

#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from textwrap import wrap
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit


# In[2]:


def make_df(path, pKdRb, cNaCl, cRbCl, cCaCl2):
    '''
    Function to make a pandas dataframe for specific variables
    input:
    output:
    df - pandas dataframe
    '''
    
    output_files_dir = os.path.join(path, pKdRb)
    
    columns = ['pHbulk',
               'fdis_A',
               'fdis_AH',
               'fdis_ANa',
               'fdis_ARb',
               'fdis_ACa',
               'fdis_A2Ca',
               'xbulk_Rb']

    if float(cRbCl) > 0:
        sysfile_dir = os.path.join(output_files_dir, 'system.cRbCl' + cRbCl + '*.dat')
    else:
        sysfile_dir = os.path.join(output_files_dir, 'system.cNaCl' + cNaCl + 'pH*.dat')
        
    # 1. Initialize an empty list to collect row dictionaries
    rows_list = []
        
    for file_name in gb.glob(sysfile_dir):
        # Initialize variables to prevent errors if a file is malformed/missing data
        pHbulk = fdis_A = fdis_AH = fdis_ANa = fdis_ACa = fdis_A2Ca = fdis_ARb = None

        with open(file_name, 'r') as input_file:
            lines = input_file.readlines()
            for line in lines:
                if 'pHbulk' in line:
                    word = line.split('=')
                    pHbulk = float(word[-1])            
                if 'avfdisA' in line:
                    word = line.split()
                    fdis_A = float(word[2])
                    fdis_AH  = float(word[3])
                    fdis_ANa = float(word[4])
                    fdis_ACa = float(word[5])
                    fdis_A2Ca = float(word[6])           
                    fdis_ARb  = float(word[7]) 
                if 'xbulk%Rb ' in line:
                    word = line.split('=')
                    xbulk_Rb = float(word[-1])

        # 2. Append the dictionary to our list instead of a DataFrame
        rows_list.append({
            'pHbulk': pHbulk,
            'fdis_A' : fdis_A,
            'fdis_AH' : fdis_AH,
            'fdis_ANa' : fdis_ANa,
            'fdis_ACa' : fdis_ACa,
            'fdis_A2Ca' : fdis_A2Ca,
            'fdis_ARb' : fdis_ARb,
            'xbulk_Rb' : xbulk_Rb 
        }) 
        
    # 3. Build the final DataFrame all at once
    df = pd.DataFrame(rows_list, columns=columns)
        
    df.sort_values('pHbulk', ascending=True, inplace=True)
    df.reset_index(drop=True, inplace=True)   # to reset index after sorting
    return df


# In[ ]:





# # analyze

# In[3]:


path_to_result = '/Users/rnap/data/data-pafiber/run_pafiber_born/18PA/'
search_dir=path_to_result+"chargeregKpK-0.6"
pKdRb_list=gb.glob(search_dir)
fig_dir="figures-chargeregKpK-0.6-test"
result_dir = os.path.join(path_to_result, 'figures', fig_dir)
os.makedirs(result_dir, exist_ok=True)


# In[4]:


cNaCl_to_analyze       = '0.0'
cNaCl_to_analyze       = '0.000' 
cCaCl2_to_analyze       = '0.0'
cRbCl_list= ['0.200','0.150','0.100','0.050','0.010','0.001']
# add a dict
data_dict={} 
for cRbCl_to_analyze  in cRbCl_list:
    data_df= make_df(path_to_result,
                     pKdRb_list[0],
                     cNaCl_to_analyze,
                     cRbCl_to_analyze,
                     cCaCl2_to_analyze)
    data_dict[cRbCl_to_analyze]=data_df


# In[5]:


# data type and size
N=180
delta=0.2
#vNa = 4.4451776756409545E-003
vK  = 1.1008442251073381E-002
vRb = 1.4710226951490483E-002
vCs= 1.9509135393555040E-002

switchRb_with_K= True
switchRb_with_Cs= False
sigmaGlu=17.3*2
pKa=5.0


# In[6]:


if switchRb_with_K :
    vRb=vK
if switchRb_with_Cs :
    vRb=vCs


# In[7]:


#data_dict['0.100']


# In[8]:


for key in data_dict.keys():
    print(key)


# In[9]:


def make_plot_fdisA(saveFig):
    fig=plt.figure()
    for key in data_dict.keys():
        data_to_plot= data_dict[key]
        pH = data_to_plot['pHbulk'].values
        fdisA = data_to_plot['fdis_A'].values
        plt.plot(pH,fdisA)
        plt.xlabel('pH',fontsize=15)
        plt.ylabel(r'$f_{A}$',fontsize=15)
    plt.show()
    if(saveFig):
        image_format = '.pdf'
        output_file = 'fdis_vs_pH' + 'cRbCl'+ image_format
        fname = os.path.join(result_dir, output_file)
        fig.savefig(fname, bbox_inches='tight', transparent=True)


# In[10]:


def make_plot_fdisA_cRbCl(key,saveFig):
    fig=plt.figure()
    data_to_plot= data_dict[key]
    pH    = data_to_plot['pHbulk'].values        
    fdisA = data_to_plot['fdis_A'].values
    plt.plot(pH,fdisA)
    plt.xlabel('pH',fontsize=15)
    plt.ylabel(r'$f_{A}$',fontsize=15)
    plt.show()
    if(saveFig):
        image_format = '.pdf'
        output_file = 'fdis_vs_pH' + 'cRbCl-pKd' + image_format
        fname = os.path.join(result_dir, output_file)
        fig.savefig(fname, bbox_inches='tight', transparent=True)


# In[11]:


def make_data_fdisA_cRbCl(key):
    data_to_plot= data_dict[key]
    pH    = data_to_plot['pHbulk'].values        
    fdisA = data_to_plot['fdis_A'].values
    data=[(x,y) for (x,y) in zip(pH,fdisA)]
    output_file = 'fdis_vs_pH' + 'cRbCl'+key+ ".dat"
    fname = os.path.join(result_dir, output_file)
    np.savetxt(fname,data)


# In[12]:


def make_ideal_fdisA(pKa,pHlist):
    """
        compute ideal dissociation 
    """
    fdisA=[]
    for pH in pHlist:
        fdis=1.0/(1.0+10**(pKa-pH))
        fdisA.append(fdis)
    return(fdisA)


# In[13]:


def save_data_fdisA_cRbCl_all():
    
    first_key = next(iter(data_dict))
    first_value = data_dict[first_key] 
    pH = data_dict[first_key]['pHbulk'].values
    fdisAideal=make_ideal_fdisA(pKa,pH)
    
    data=[]
    data.append(pH)
    data.append(fdisAideal)
    
    header_txt="pH, fdis_ideal"
    
    for key in data_dict.keys(): # loop salt 
        data_to_save= data_dict[key]
        data.append(data_to_save["fdis_A"].values)
        header_txt=header_txt+", c="+str(key)+"M"
    
    # zip(*fdis_state) unpacks the rows and groups them column-wise
    data_out = [list(column) for column in zip(*data)]
    #print(header_txt)
    
    output_file = 'fdisA_vs_pH_cRbCl'+'.txt'
    fname = os.path.join(result_dir, output_file)
    
    np.savetxt(fname, data_out, header=header_txt)
    


# In[14]:


def make_plot_disociat_cRbCl(key):
    states = {"A" : 0, "AH" : 1, "ANa" : 2, "ARb" : 3,"ACa": 4, "A2Ca":5}
    data_to_plot=data_dict[key]
    pH     = data_to_plot['pHbulk'].values
    fdis=[]
    for state in states:
        fdis.append(data_to_plot["fdis_"+state].values)
        plt.plot(pH,fdis[states[state]],label=state)
    plt.xlabel('pH',fontsize=15)
    plt.ylabel(r'$f_{k}$',fontsize=15)
    plt.legend()    


# In[15]:


def make_data_disociat_cRbCl(key):
    states = {"A" : 0, "AH" : 1, "ANa" : 2, "ARb" : 3,"ACa": 4, "A2Ca":5}
    data_to_plot=data_dict[key]
    pH     = data_to_plot['pHbulk'].values
    fdis=[]
    y=[]
    for state in states:
        fdis_state=data_to_plot["fdis_"+state].values
        data=[(x,y) for (x,y) in zip(pH,fdis_state)]
        output_file = 'fdis_'+state+'_vs_pH_cRbCl'+key+".txt"
        fname = os.path.join(result_dir, output_file)
        np.savetxt(fname,data)


# In[16]:


key="0.100"
#make_plot_fdisA_cNaCl(key,False)
#make_data_fdisA_cNaCl(key)
for key in data_dict.keys():
    make_data_fdisA_cRbCl(key)


# In[17]:


save_data_fdisA_cRbCl_all() 


# In[18]:


saveFig=False
make_plot_fdisA(saveFig)


# In[19]:


key='0.100'
make_plot_disociat_cRbCl(key)
make_data_disociat_cRbCl(key)


# In[20]:


for key in data_dict.keys():
    print(key)


# In[21]:


def make_data_ionexcess_Rb_condensed(key):
    """
        Calculate for dataframe data_dict and key = salt concentration the  
        total number of ion per unit length that are condensed
        for C16-V2A2E2 with line denisty 17.3 1/nm == 18PA 
    """
    # key = salt concentration
    states = {"A" : 0, "AH" : 1, "ANa" : 2, "ARb" : 3,"ACa": 4, "A2Ca":5}
    data=data_dict[key]
    # get pH 
    pH  = data['pHbulk'].values
    fdis=[]
    state="ARb"
    fdis=data["fdis_"+state].values
    gamma_cond=fdis*17.3*2 # explicit for 17.3 GLu per nm and 2 
    # total number of ion per uni lenght that are condensed
    
    # output 
    return(pH,gamma_cond)


# In[22]:


def make_delta_A(delta,N,R=0):
    """
        Compute geometrical factor for cylinder symmetry
    """
    deltaA=np.zeros(N)
    for i in range(N):
        rmid= (i+0.5) * delta + R # in fortran rmid= (i-0.5) * delta + R , becuase start i=1 here i=0
        deltaA[i]= 2.0 * np.pi * rmid * delta
    return(deltaA)       


# In[23]:


def read_density_Glu(path,N):
    """
        reads GLu density in from file 
    """
    file_name=path+"/rhoEpa.in"
    rho_Glu=np.zeros(N)
    with open(file_name, 'r') as input_file:
        lines = input_file.readlines()
        count=0
        for line in lines:
            word = line.split()
            rho = float(word[1])
            if count < N:
                rho_Glu[count]=rho
            count=count+1    
    return(rho_Glu) 


# In[24]:


def read_iondensity_Rb(path,cRbCl,pH,volRb):
    """
        reads volume fraction of Rb ion 
        input: 
            path : folder location
            cRbCl : concentration cRbCl
            pH   : pH
            volRb : volume Rb ion
        return: 
            rhoRb : density Rb
            
    """
    file_name=path+"/xRbions.cRbCl"+str(cRbCl)+"pH"+str(f"{pH:.3f}")+".dat"
    data = np.loadtxt(file_name)
    rhoRb = data[:, 1]/volRb # Second column
    return(rhoRb)


# In[25]:


def make_data_ionexcess_Rb(key,delta,N,volRb):
    """
        Calculates for the conditon of dataframe data_dict with key = salt concentration the  
        total number of excess free ions per unit length 
        for C16-V2A2E2 with line denisty 17.3 1/nm == 18PA 
        input: 
            key =salt concentration
            delta = lattice spacing
            N = lenght iondensities
            volRb =volume Rb
        return:
            pHval,gamma_condensed, gamma_free, gamma_tot
            
    """
    
    # key = salt concentration
    data=data_dict[key]
    # get pH 
    pHval  = data['pHbulk'].values
    # xbulk 
    xbulkRbval  = data['xbulk_Rb'].values
    
    
    # total number of ion per unit length that are condensed
    state="ARb"
    fdis=data["fdis_"+state].values
    gamma_condensed=fdis*17.3*2
   
    # ion excess free ion  
    
    gamma_free=[]
    gamma_tot=[]
    
    # 0. need deltaA geometric factor 
    deltaA = make_delta_A(delta,N)
    
    # 1. read Glu density 
    rho_Glu = read_density_Glu(pKdRb_list[0],N)
    
    
    for pH,xbulk  in zip(pHval,xbulkRbval):
      
        # 2. read in Rb+ ion density
    
        rhoRb = read_iondensity_Rb(pKdRb_list[0],key,pH,volRb)
    
        # 3. preform integral \int A(r) \rho_{Na^+}= sum_ delta A(i) \rho_{N^+}(i)
        
        gamma=0.0
        ref=0.0
        nGlu =0.0
        for i in range(len(rhoRb)):
            gamma = gamma + deltaA[i] * rhoRb[i]
            ref = ref + deltaA[i] *  xbulk/volRb
            nGlu= nGlu + deltaA[i] * rho_Glu
               
        gamma_free.append(gamma-ref)
    
    gamma_tot=gamma_free+gamma_condensed
                       
    # output
    return pHval,gamma_condensed, gamma_free, gamma_tot    
    


# In[26]:


def make_plot_ionexcess_Rb(key,saveFig):
   
    pH, gamma_cond, gamma_free, gamma_tot= make_data_ionexcess_Rb(key,delta,N,vRb)

    fig=plt.figure()
    plt.plot(pH,gamma_cond,label='condensed')
    plt.plot(pH,gamma_free,label="free")
    plt.plot(pH,gamma_tot,label="total")    
    plt.xlabel('pH',fontsize=15)
    plt.ylabel(r'$\Gamma$',fontsize=15)
    plt.legend() 
    plt.show()
    
    if(saveFig):
        image_format = '.pdf'
        output_file = 'gamma_vs_pH' + 'cRbCl'+key + image_format
        fname = os.path.join(result_dir, output_file)
        fig.savefig(fname, bbox_inches='tight', transparent=True)
    


# In[27]:


def make_plot_ionexcess_Rb_all(saveFig):

    fig=plt.figure()
    for key in data_dict.keys():
        pH, gamma_cond, gamma_free, gamma_tot= make_data_ionexcess_Rb(key,delta,N,vRb)
        plt.plot(pH,gamma_tot,label="c="+str(key))    
    plt.xlabel('pH',fontsize=15)
    plt.ylabel(r'$\Gamma$',fontsize=15)
    plt.legend() 
    plt.show()
    
    if(saveFig):
        image_format = '.pdf'
        output_file = 'gamma_vs_pH' + 'cRbCl' + image_format
        fname = os.path.join(result_dir, output_file)
        fig.savefig(fname, bbox_inches='tight', transparent=True)
    


# In[28]:


def save_data_ionexcess_Rb(key):
   
    pH,gamma_cond, gamma_free, gamma_tot = make_data_ionexcess_Rb(key,delta,N,vRb)
    data=[(x,y,z,u) for (x,y,z,u) in zip(pH,gamma_cond, gamma_free, gamma_tot)]
          
    output_file = 'gamma_vs_pH' + 'cRbCl'+key +'.txt'
    fname = os.path.join(result_dir, output_file)
    
    np.savetxt(fname, data, header='pH,gamma_cond, gamma_free, gamma_tot')
    


# In[29]:


def save_data_ionexcess_Rb_all():
    
    first_key = next(iter(data_dict))
    first_value = data_dict[first_key] 
    pH = data_dict[first_key]['pHbulk'].values
    data=[]
    data.append(pH)
    header_txt="pH, "
    
    for key in data_dict.keys():
        pH,gamma_cond, gamma_free, gamma_tot = make_data_ionexcess_Rb(key,delta,N,vRb)
        data.append(gamma_tot)
        header_txt=header_txt+", c="+str(key)+"M"
    
    # zip(*fdis_state) unpacks the rows and groups them column-wise
    data_out = [list(column) for column in zip(*data)]

    output_file = 'gamma_vs_pH' + 'cRbCl.txt'
    fname = os.path.join(result_dir, output_file)
    
    np.savetxt(fname, data_out, header=header_txt)


# In[32]:


key='0.050'
saveFig= False
make_plot_ionexcess_Rb(key,saveFig)
save_data_ionexcess_Rb(key)
key='0.100'
save_data_ionexcess_Rb(key)
make_plot_ionexcess_Rb_all(saveFig)
save_data_ionexcess_Rb_all()


# In[ ]:





# In[ ]:




