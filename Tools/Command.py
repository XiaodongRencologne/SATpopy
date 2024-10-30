import numpy as np

#%%
class get_current():
    def __init__(self,
                 object,
                 source,
                 accuracy= -100,
                 auto_convergence=True,
                 convergence_on_scatterer = None,
                 convergence_on_output_grid = None,
                 convergence_on_expansion_grid = None,
                 max_bisections = 5,
                 integration_grid_limit = 'off',
                 obsolete_conv_on_po_grid = None
                 ):
        ''' source is a python list '''
        self.object = object
        self.source = source
        self.accuracy = accuracy
        self.auto_convergence =auto_convergence
        self.convergence_on_scatterer = convergence_on_scatterer
        self.convergence_on_output_grid = convergence_on_output_grid
        self.convergence_on_expansion_grid =convergence_on_expansion_grid 
        self.max_bisections = max_bisections
        self.integration_grid_limit = integration_grid_limit
        self.obsolete_conv_on_po_grid = obsolete_conv_on_po_grid

        self.Str = self._get_str()

    def _get_str(self):
        Str = ''
        Str += 'COMMAND OBJECT ' + self.object.name + ' get_currents &\n'
        Str += '(   source                : sequence('
        for item in self.source:
            Str += 'ref('+item.name+'),'
        Str += '),&\n'
        Str += 'field_accuracy           : ' + str(self.accuracy) +',&\n'
        if self.auto_convergence:
            Str += 'auto_convergence_of_po    : ' + 'on' +',&\n'
            if self.convergence_on_scatterer != None:
                Str += 'convergence_on_scatterer  : ' + 'sequence('
                for item in self.convergence_on_scatterer:
                    Str += 'ref('+item.name +'),'
                Str += '),&\n'
            if self.convergence_on_output_grid != None:
                Str += 'convergence_on_output_grid  : ' + 'sequence('
                for item in self.convergence_on_output_grid:
                    Str += 'ref('+item.name +'),'
                Str += '),&\n'
            if self.convergence_on_expansion_grid != None:
                Str += 'convergence_on_expansion_grid  : ' + 'sequence('
                for item in self.convergence_on_expansion_grid:
                    Str += 'ref('+item.name +'),'
                Str += '),&\n'
        else:
            Str += 'auto_convergence_of_po   : ' + 'off' +',&\n'
        
        Str += 'max_bisections            : ' + str(self.max_bisections) +',&\n'
        Str += 'integration_grid_limit    : ' + self.integration_grid_limit +',&\n'
        
        if self.obsolete_conv_on_po_grid != None:
            Str += 'obsolete_conv_on_po_grid   : '
        Str += ')&\n\n'
        return Str

# %%
class get_field():
    def __init__(self,
                 object,
                 source,
                 ):
        ''' source is a python list '''
        self.object = object
        self.source = source
        self.Str = self._get_str()

    def _get_str(self):
        Str = ''
        Str += 'COMMAND OBJECT ' + self.object.name + ' get_field &\n'
        Str += '(  source                : sequence('
        for item in self.source:
            Str += 'ref('+ item.name +'),'
        Str += ')&\n'
        Str += ') &\n\n'
        return Str