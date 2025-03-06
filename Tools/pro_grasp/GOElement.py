import numpy as np
#import W2GRASP_GO, W2GRASP_EO
from .W2GRASP_GO import write_coord,write_simple_lens,write_aperture,write_scatter_cluster
from .W2GRASP_GO import write_rim,write_BoR_Mesh
from .GaussianOptics import ThinLens


'''
Define Geometrical objectors, coordinate system, simple lens.
'''
class _global_coord_sys():
    def __init__(self):
        self.origin=np.zeros((3,1))
        self.origin_g=np.zeros((3,1))
        self.name = 'mechanical_axis'
        self.Str = write_coord('mechanical_axis',[0,0,0],[0,0,0])
    
global_coord=_global_coord_sys()

class coor_sys():
    def __init__(self,origin, angle, ref_coor = global_coord, name='coor_sys'):
        self.name = name
        self.ref_coor = ref_coor
        #self.ref_coor = ref_coor
        self.origin = origin
        self.angle = angle

        self.Str = write_coord(self.name,
                               self.origin,
                               self.angle,
                               base = self.ref_coor.name)
        self.List = {'name': 'coor_feed_ref',
                     'origin':[0,0,0],
                     'angle': [0,0,0],
                     'base':'mechanical_axis.cs',
                     'Type':'coor_sys'}
        
class reflector():
    def __init__(self,):
        pass
class simple_lens():
    def __init__(self, coor_sys, 
                 diameter, 
                 refractive_index, loss_tangent = 0, 
                 r1 = None, r2 = None, bs1 = 0, bs2 = 0,
                 thickness = '0',
                 surf_f1 = '',
                 surf_f2 = '',
                 lengthUnit = 'mm',
                 coating_surf1 = '', coating_surf2 = '',
                 name='simple_lens'):
        self.name = name
        self.coor_name = coor_sys.name
        self.diameter = diameter
        self.refractive_index = refractive_index
        self.loss_tangent = loss_tangent
        self.r1 =r1
        self.r2 = r2
        self.bs1 = bs1
        self.bs2 = bs2
        self.thickness =thickness
        self.surf_f1 = surf_f1
        self.surf_f2 = surf_f2
        self.lengthUnit = lengthUnit
        self.coating_surf1 = coating_surf1
        self.coating_surf2 = coating_surf2

        self.Str = write_simple_lens(self.name,
                                     self.coor_name,
                                     str(self.diameter) + ' cm',
                                     self.refractive_index,
                                     self.loss_tangent,
                                     r1 = self.r1, r2 = self.r2, 
                                     bs1 = self.bs1, bs2 = self.bs2,
                                     thickness = self.thickness,
                                     surf1_file = self.surf_f1,surf2_file = self.surf_f2,
                                     lengthUnit_file = self.lengthUnit,
                                     coating_surf1 = self.coating_surf1, coating_surf2 = self.coating_surf2)
        

class Aperture_screen():
    def __init__(self, coor_sys, 
                 rim,
                 infinity_shadow = 'on',
                 name='aperture'):
        self.name = name
        self.coor_name = coor_sys.name
        self.rim = rim.name
        self.infinity_shadow = infinity_shadow
        self.Str = write_aperture(self.name,self.coor_name,self.rim, self.infinity_shadow)

class rim():
    def __init__(self, center, 
                 side_lengths, 
                 Type = 'rectangular_rim',name='rim'):
        '''Type = 'rectangular_rim, elliptical_rim'''
        self.center = center
        self.side_lengths = side_lengths
        self.Type = Type
        self.name = name
        self.Str = write_rim(self.name,self.center,self.side_lengths, Type = self.Type)


class scatterer_cluster():
    def __init__(self,scatter_list,name = 'cluster'):
        self.name = name
        self.scatters = scatter_list
        self.Str = write_scatter_cluster(self.scatters)

class BoR_MoM():
    def __init__(self,
                 coor_sys,
                 region,
                 nodes,
                 linear_segments,
                 #cubic_segments,
                 length_unit = 'mm',
                 advanced_regions='',
                 name='BoR_Mesh'):
        self.name =name
        self.coor_name = coor_sys.name
        self.region = region
        self.nodes = nodes
        self.linear_segments = linear_segments
        self.length_unit = length_unit
        self.Str = write_BoR_Mesh(name,
                                  self.coor_name,
                                  self.region,
                                  self.nodes,
                                  self.linear_segments,
                                  #cubic_segments,
                                  self.length_unit,
                                  advanced_regions=advanced_regions)
        

class BoR_MoM_lens():
    def __init__(self,
                 coor_sys,
                 lens_index,
                 lens_t,
                 layer_index,
                 AR_t,
                 surface1,
                 surface2,
                 length_unit = 'mm',
                 advanced_regions='',
                 name='BoR_Mesh'):
        self.name =name
        self.coor_name = coor_sys.name
        self.lens_index= lens_index
        self.layer_index = layer_index
        region = [[1,lens_index**2,1,0]]
        n =2
        for item in layer_index:
            region.append([n, item**2 ,1,0])
            n+=1
        self.region = region

        # define nodes
        N_surface = len(layer_index) +2
        rho = np.repeat(surface1[:,0],N_surface).reshape(-1,N_surface).T.ravel()
        #print(rho,rho.size)
        surf = np.concatenate((surface1[:,1],-surface2[:,1]+lens_t))

        L_t1 = 0
        L_t2 = lens_t+0
        for n in range(len(AR_t)):
            if n%2 ==0:
                L_t1 += AR_t[n]
                surf = np.concatenate((surf,surface1[:,1]-L_t1))
            else:
                L_t2 += AR_t[n]
                surf = np.concatenate((surf,-surface2[:,1]+L_t2))
        N_nodes = np.array(list(range(1,rho.size+1)))
        #print(N_nodes.size,rho.size,surf.size)
        self.nodes = np.concatenate((N_nodes,rho, surf)).reshape(3,-1).T

        # define order of nodes
        node_A = np.array([])
        node_B = np.array([])
        region1 = np.array([])
        region2 = np.array([])
        N_nodes_list = []
        for n in range(2*len(layer_index)+2):
            if n%2 ==0:
                N_nodes_list.append(int(surface1.size/2))
            else:
                N_nodes_list.append(int(surface2.size/2))
        n_start = 0
        End = []
        for n in range(N_surface):
            if n%2 ==0:
                node_A = np.append(node_A, N_nodes[n_start:n_start+N_nodes_list[n]-1])
                node_B = np.append(node_B, N_nodes[n_start+1:n_start+N_nodes_list[n]])
                if n == 0:
                    region1 = np.append(region1,np.zeros(N_nodes_list[n]-1)+(n+1))
                    region2 = np.append(region2,np.zeros(N_nodes_list[n]-1)+(n+2))
                elif n == N_surface-2:
                    region1 = np.append(region1,np.zeros(N_nodes_list[n]-1)+(n))
                    region2 = np.append(region2,np.zeros(N_nodes_list[n]-1))
                else:
                    region1 = np.append(region1,np.zeros(N_nodes_list[n]-1)+(n))
                    region2 = np.append(region2,np.zeros(N_nodes_list[n]-1)+(n+2))
                n_start +=N_nodes_list[n]
            else:
                node_A = np.append(node_A, N_nodes[n_start:n_start+N_nodes_list[n]-1])
                node_B = np.append(node_B, N_nodes[n_start+1:n_start+N_nodes_list[n]])
                if n == N_surface-1:
                    region1 = np.append(region1,np.zeros(N_nodes_list[n]-1))
                else:
                    region1 = np.append(region1,np.zeros(N_nodes_list[n]-1)+(n+2))
                region2 = np.append(region2,np.zeros(N_nodes_list[n]-1)+(n))
                n_start +=N_nodes_list[n]
            End.append(node_B[-1])
        # lens
        node_A = np.append(node_A, End[0])
        node_B = np.append(node_B, End[1])   
        region1 = np.append(region1,1)
        region2 = np.append(region2,0) 
        for n in range(len(layer_index)):
            if n%2 == 0:
                node_A = np.append(node_A, End[n+2])
                node_B = np.append(node_B, End[n])
                region1 = np.append(region1,n+2)
            else:
                node_A = np.append(node_A, End[n])
                node_B = np.append(node_B, End[n+2])
                region1 = np.append(region1,n+2)
            region2 = np.append(region2, 0)
        self.linear_segments = np.concatenate((np.array(list(range(1,node_A.size+1))),
                                               node_A,node_B,region1, region2,-np.ones(node_A.size),np.zeros(node_A.size))).reshape(7,-1).T
        self.length_unit = length_unit
        self.Str = write_BoR_Mesh(name,
                                  self.coor_name,
                                  self.region,
                                  self.nodes.tolist(),
                                  self.linear_segments.tolist(),
                                  #cubic_segments,
                                  self.length_unit,
                                  advanced_regions=advanced_regions)



# Quasi-optics reflector, off-axis, conic mirror or polynomial mirror. 

class Quasi_optics_reflector():
    def __init__(self,
                 win,din,f,Angle,Lambda,
                 align_type='BW',
                 align_coord='',
                 name='off_axis_reflector'):
        self.Gaussian_params = Thinlens(win,din,f,Lambda)
        self.r1 = Mirror['Rin']
        self.r2 = Mirror['Rout']
        self.dout = Mirror['dout']
        self.din = din
        self.angle= Angle
        
        pass