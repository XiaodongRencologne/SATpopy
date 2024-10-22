#%%
import numpy as np
from numpy.polynomial.polynomial import polyval
import matplotlib.pyplot as plt
# %%
def EvenAsphere(R,k,coeffi_even):
    c =1/R
    coeff = np.append(np.array([0.0]),np.array(coeffi_even))
    def surf1(rho):
        z = c*rho**2/(1+np.sqrt(1-(1+k)*c**2*rho**2))   
        rho2 =rho**2
        z += polyval(rho2,coeff)
        return z
    return surf1
def zemax2RSF(Np,Kspace,Ktip,lens_para):
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
    with open(lens_para['name']+'.rsf','w') as f:
        f.writelines(lens_para['name']+'\n')
        f.writelines(str(Np)+', '+str(Kspace)+', '+str(Ktip)+'\n')
    D = lens_para['D']
    if lens_para['type'] == 'EvenAsphere':
        surf_fuc = EvenAsphere(lens_para['R'],lens_para['K'],
                               lens_para['co'])
    if Kspace == 0:
        rho = np.linspace(0,D/2,Np)
        z = surf_fuc(rho)
        with open(lens_para['name'] + '.rsf','a') as f:
            f.writelines(str(rho.min())+', '+str(rho.max()) +'\n')
            for n in range(Np):
                f.writelines(str(z[n])+'\n')
# %%
f1= EvenAsphere(39.996,1.891,[-0.013,0,0])
f2= EvenAsphere(-41.313,-0.730,[0.013,1.379*10**(-6),1.83*10**(-8)])

# %%
lens3_face1 = {'R':3.999622611956000E+001,
               'K':1.890831602573785E+000,
               'type':'EvenAsphere',
               'co': [-1.280526111650124E-002, 0, 0],
               'D': 1.979846446331175E+001*2,
               'name':'lens3_f1'}
lens3_face2 = {'R':-4.131265169623969E+001,
               'K':-7.301384185400588E-001,
               'type':'EvenAsphere',
               'co': [1.322423720533787E-002,1.379457394858853E-006,1.829653402510510E-008],
               'D': 1.979846446331175E+001*2,
               'name':'lens3_f2'}

lens2_face1 = {'R':39.999,
               'K':-7.465,
               'type':'EvenAsphere',
               'co': [-3.324674200128043E-003, 4.339567195952648E-006, -6.609342622830201E-009],
               'D': 2.221349434353119E+001*2,
               'name':'lens2_f1'}
lens2_face2 = {'R':39.999,
               'K':-4.733,
               'type':'EvenAsphere',
               'co': [-7.252526369706933E-003,0,0],
               'D': 2.221349434353119E+001*2,
               'name':'lens2_f2'}

lens1_face1 = {'R':7.317849839565703E+001,
               'K':-7.465,
               'type':'EvenAsphere',
               'co': [-8.194722933835832E-004, 5.647849880421205E-006, -8.282440810597639E-009],
               'D': 2.220424646491709E+001*2,
               'name':'lens1_f1'}
lens1_face2 = {'R':39.999,
               'K':-4.733,
               'type':'EvenAsphere',
               'co': [4.444393059968795E-003, 2.091776374381418E-006, -6.822315879457812E-009],
               'D': 2.220424646491709E+001*2,
               'name':'lens1_f2'}
# %%
lenses = [lens1_face1,lens1_face2,
          lens2_face1,lens2_face2,
          lens3_face1,lens3_face2]
# %%
Np = 201
Kspace = 0
Ktip = 0
for item in lenses:
    zemax2RSF(Np,Kspace,Ktip,item)
# %%
