#%%
import numpy as np
from numpy.polynomial.polynomial import polyval
import matplotlib.pyplot as plt
# %%
def EvenAsphere(R,k,coeffi_even):
    c =1/R
    coeff = np.append(np.array([0.0]),np.array(coeffi_even))
    print(coeff)
    def surf1(rho):
        z = c*rho**2/(1+np.sqrt(1-(1+k)*c**2*rho**2))   
        rho2 =rho**2
        z += polyval(rho2,coeff)
        return z
    return surf1
def zemax2RSF(Np,Kspace,Ktip,lens_para,outputfolder='',sign = 1):
    '''
    Rotationally symmetric surface.
    It is one-demensional surface which is a function of radial Rho. 
    rho =sqrt(x^2+y^2)

    Lens_para = {'R': 500,
                 'K': -2.1,
                 'type': 'EvenAsphere',
                 'co': [1,2,3],
                 'D' : 200,
                 'name':'lens_face1'}
    '''
    with open(outputfolder+lens_para['name']+'.rsf','w') as f:
        f.writelines(lens_para['name']+'\n')
        f.writelines(str(Np)+' '+str(Kspace)+' '+str(Ktip)+'\n')
    D = lens_para['r']*2
    if lens_para['type'] == 'EvenAsphere':
        surf_fuc = EvenAsphere(lens_para['R'],lens_para['K'],
                               lens_para['co'])
    if Kspace == 1:
        rho = np.linspace(0,D/2,Np)
        z = sign * surf_fuc(rho)
        data = np.append(rho,z).reshape(2,-1).T
        with open(outputfolder+lens_para['name'] + '.rsf','a') as f:
            #f.writelines(str(rho.min())+' '+str(rho.max()) +'\n')
            np.savetxt(f,data,delimiter=' ',fmt = '%10.8f')
            """
            for n in range(Np):
                f.writelines(str(rho[n]) + ' ' +str(z[n])+'\n')
            """
    return rho, z

# %%
Diameter = 22.4
lens1_face1 = {'R':1.035545854703427E+002,
               'K':-3.000122786911282E+001,
               'type':'EvenAsphere',
               'co': [1.450211740738392E-003, 2.848581234586182E-006, -5.021767369327404E-009],
               'r': Diameter,
               'name':'lens1_f1'}
lens1_face2 = {'R':-1.672578325120606E+002,
               'K':2.999835005933332E+001,
               'type':'EvenAsphere',
               'co': [6.860967807195721E-003, 4.447195339970823E-007, -2.621780005775388E-009],
               'r': Diameter,
               'name':'lens1_f2'}
lens2_face1 = {'R':8.189215891762041E+001,
               'K':-3.000046666787026E+001,
               'type':'EvenAsphere',
               'co': [1.564605534523395E-003, 3.063497794251431E-006, -5.002469546126761E-009],
               'r': Diameter,
               'name':'lens2_f1'}
lens2_face2 = {'R':4.173343843027721E+001,
               'K':-1.046030791067757E+000,
               'type':'EvenAsphere',
               'co': [-8.382357913364686E-003,-1.390480360528942E-006,1.846691641792131E-010],
               'r': Diameter,
               'name':'lens2_f2'}

lens3_face1 = {'R':-4.023494153544174E+001,
               'K':-3.000640715654068E+001,
               'type':'EvenAsphere',
               'co': [4.058867092144052E-004, 5.000810092877716E-006, 5.006698115075333E-009],
               'r': Diameter,
               'name':'lens3_f1'}
lens3_face2 = {'R':-3.997787293605576E+001,
               'K':-1.818966462620520E+001,
               'type':'EvenAsphere',
               'co': [2.290398386941478E-003,5.008243629130054E-006,1.945114042732393E-009],
               'r': Diameter,
               'name':'lens3_f2'}
# %%
lenses = [lens1_face1,lens1_face2,
          lens2_face1,lens2_face2,
          lens3_face1,lens3_face2]
lens_inf = np.cumsum(np.array([0.0, 3.5,51.524,3.756,10.008,2.009]))
# %%
Np = 101
Kspace = 1
Ktip = 0
n=0
fig = plt.figure(figsize=(10,10))

outputfolder = 'output/srf2/'

for item in lenses:
    print(n)
    print(item)
    if n%2 ==0:
        x, z = zemax2RSF(Np,Kspace,Ktip,item,outputfolder=outputfolder,sign = 1)
        plt.plot(x, z+lens_inf[n],'k-')
        plt.plot(-x, z+lens_inf[n],'k-')
    else:
        x, z = zemax2RSF(Np,Kspace,Ktip,item,outputfolder=outputfolder,sign = -1)
        plt.plot(x, -z+lens_inf[n],'k-')
        plt.plot(-x, -z+lens_inf[n],'k-')
    
    n+=1
plt.axis('equal')
plt.xlim([-25,25])
plt.show()

# %%
