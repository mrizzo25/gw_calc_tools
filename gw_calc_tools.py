
# coding: utf-8

# In[ ]:

'''
Web tool for calculating salient quantities in
the analysis of astrophysical systems that produce
gravitational waves.
'''


# In[32]:

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import lalsimulation as lalsim
get_ipython().magic(u'matplotlib inline')

MSUN_SI = 1.98847e30
MRSUN_SI=.476625061404649406193430731479084713e3   #MSUN in geometrized units
C_SI = 2.98e8
G_SI = 6.674e-11


# In[2]:

#Define dictionary to hold parameters and values
params={'m1':1.5,'m2':1.5,'S1':[0.5,.5,.5],'S2':[.5,.5,.5],'L':[.5,.5,.5],'eosname':'MPA1'}


# In[46]:

#Returns mass ratio
def q(params):
    #Unpack params
    m1=params['m1']
    m2=params['m2']
    
    q=m2/m1
    return q


# In[5]:

#Returns chirp mass given component masses
def chirpMass(params):
    #Unpack params
    m1=params['m1']
    m2=params['m2']
    
    Mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    return Mc


# In[6]:

#Returns symmetric mass ratio given component masses
def eta(params):
    #Unpack params
    m1=params['m1']
    m2=params['m2']
    
    eta = m1*m2/(m1+m2)
    return eta


# In[56]:

#Returns mass-weighted spin given spins, component masses
#Need to check units on these ones
def xiEff(params):
    #Unpack params
    m1=params['m1']
    m2=params['m2']
    S1=params['S1']
    S2=params['S2']
    L=params['L']
    
    xiEff=0
    M=m1+m2
    Lmag=(L[0]**2.+L[1]**2.+L[2]**2.)**(1./2.)
    Lhat=[L[0]/Lmag,L[1]/Lmag,L[2]/Lmag]
    for i in range(len(L)):
        xiEff+=(S1[i]/m1)*Lhat[i] + (S2[i]/m2)*Lhat[i]
    
    return (1./M)*xiEff


# In[57]:

#Returns effective precession spin parameter
def xiP(params):
    #Unpack params
    m1=params['m1']
    m2=params['m2']
    S1=params['S1']
    S2=params['S2']
    L=params['L']
    
    q=m2/m1
    Lmag=(L[0]**2.+L[1]**2.+L[2]**2.)**(1./2.)
    Lhat=[L[0]/Lmag,L[1]/Lmag,L[2]/Lmag]
    
    #Calculate component of spins in direction of orbital plane
    S1perp=[S1[1]*Lhat[2]-S1[2]*Lhat[1],-(S1[0]*Lhat[2]-S1[2]*Lhat[0]),S1[0]*Lhat[1]-S1[1]*Lhat[0]]
    S2perp=[S2[1]*Lhat[2]-S2[2]*Lhat[1],-(S2[0]*Lhat[2]-S2[2]*Lhat[0]),S2[0]*Lhat[1]-S2[1]*Lhat[0]]
    S1perpMag=(S1perp[0]**2. + S1perp[1]**2. + S1perp[2]**2.)**(1./2.)
    S2perpMag=(S2perp[0]**2. + S2perp[1]**2. + S2perp[2]**2.)**(1./2.)
    B1=2.+3.*q/2.
    B2=2.+3./(2.*q)
    
    xiP=1./(B1*m1**2.) * max(B1*S1perpMag,B2*S2perpMag)
    
    return xiP


# In[41]:

#Returns tidal deformability parameter given mass and EOS
def Lambdas(params):
    #Unpack params
    eosname=params['eosname']
    m1=params['m1']
    m2=params['m2']
    
    #Create EOS/Family structures
    eos=lalsim.SimNeutronStarEOSByName(eosname)
    fam=lalsim.CreateSimNeutronStarFamily(eos)
    r1=lalsim.SimNeutronStarRadius(m1*MSUN_SI,fam)
    r2=lalsim.SimNeutronStarRadius(m2*MSUN_SI,fam)
    k1=lalsim.SimNeutronStarLoveNumberK2(m1*MSUN_SI,fam)
    k2=lalsim.SimNeutronStarLoveNumberK2(m2*MSUN_SI,fam)
    c1=m1*MRSUN_SI/r1
    c2=m2*MRSUN_SI/r2
    Lambda1= (2.0/3.0) * k1 / c1**5.
    Lambda2= (2.0/3.0) * k2 / c2**5.
    
    return [Lambda1,Lambda2]


# In[42]:

Lambdas(params)


# In[ ]:



