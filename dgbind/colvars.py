from jinja2 import Template
import MDAnalysis as mda
import os
import StringIO

class ColvarRMSD:
    _template = """
colvar {
    name {{name}}

    rmsd {
        atoms {
            atomnumbers {
                {%- for column in selected_atoms|batch(10) %}
                {% for atom in column -%}
                {{ "%6i" | format(atom.number+1) }} 
                {%- endfor -%}
                {% endfor %}
            }
        }
        refPositionsFile {{ inputdir }}/{{ refpdb }}
        refPositionsCol  B
    }
}

harmonic {
    name          {{ name }}
    colvars       {{ name }}
    forceconstant {{ forceconstant|default(10.0) }}
    centers       {{ refvalue|default(0.0) }}
}
"""

    def __init__(self, name, universe, **kwargs):
        self.template = Template(self._template)
        self.spec = kwargs
        self.spec['name'] = name
        if self.spec.has_key('refpdb'):
            universe.load_new(self.spec['refpdb'])
        if self.spec.has_key('receptor'): self.spec['selection'] = self.spec['selection'].replace('receptor', self.spec['receptor'])
        if self.spec.has_key('ligand'): self.spec['selection'] = self.spec['selection'].replace('ligand', self.spec['ligand'])
        self.spec['selected_atoms'] = universe.select_atoms(self.spec['selection'])
        self.spec['selected_atoms'].set_bfactors(1)
        self.spec['refpdb'] = '%s.ref' % name

    def write(self):
        return self.template.render(self.spec)

class ColvarAngle:
    _template = """
colvar {
    name {{name}}

    {{ angletype }} {
        {%- for selected_atoms in refatoms %}
        group{{ loop.index }} {
            atomnumbers {
                {%- for column in selected_atoms|batch(10) %}
                {% for atom in column -%}
                {{ "%6i" | format(atom.number+1) }} 
                {%- endfor -%}
                {% endfor %}
            }
        }
        {%- endfor %}
    }
}

harmonic {
    name          {{ name }}
    colvars       {{ name }}
    forceconstant {{ forceconstant|default(0.1) }}
    centers       {{ centers }}
}
"""

    def __init__(self, name, universe, **kwargs):
        self.template = Template(self._template)
        self.spec = kwargs
        self.spec['name'] = name
        self.spec['refatoms'] = [universe.select_atoms(self.spec['refatoms'][ref]) for ref in self.spec['angle']]
        self.spec['angletype'] = 'angle' if len(self.spec['refatoms']) == 3 else 'dihedral'

    def write(self):
        return self.template.render(self.spec)

class ColvarDistance:
    _template = """
colvar {
    name {{name}}

    distance {
        {%- for selected_atoms in refatoms %}
        group{{ loop.index }} {
            atomnumbers {
                {%- for column in selected_atoms|batch(10) %}
                {% for atom in column -%}
                {{ "%6i" | format(atom.number+1) }} 
                {%- endfor -%}
                {% endfor %}
            }
        }
        {%- endfor %}
    }
}

harmonic {
    name          {{ name }}
    colvars       {{ name }}
    forceconstant {{ forceconstant|default(10.0) }}
    centers       {{ centers }}
}
"""

    def __init__(self, name, universe, **kwargs):
        self.template = Template(self._template)
        self.spec = kwargs
        self.spec['name'] = name
        self.spec['refatoms'] = [universe.select_atoms(self.spec['refatoms'][ref]) for ref in self.spec['distance']]

    def write(self):
        return self.template.render(self.spec)

class ColvarOmega:
    _template = """
colvar {
    name {{name}}

    orientation {
        atoms {
            atomnumbers {
                {%- for column in selected_atoms|batch(10) %}
                {% for atom in column -%}
                {{ "%6i" | format(atom.number+1) }} 
                {%- endfor -%}
                {% endfor %}
            }
        }
        refPositionsFile {{ inputdir }}/{{ refpdb }}
        refPositionsCol  B
    }
}

harmonic {
    name          {{ name }}
    colvars       {{ name }}
    forceconstant {{ forceconstant|default(500) }}
    centers       {{ centers|default("(1.0, 0.0, 0.0, 0.0)") }}
}
"""

    def __init__(self, name, universe, **kwargs):
        self.template = Template(self._template)
        self.spec = kwargs
        self.spec['name'] = name
        if self.spec.has_key('receptor'): self.spec['selection'] = self.spec['selection'].replace('receptor', self.spec['receptor'])
        if self.spec.has_key('ligand'): self.spec['selection'] = self.spec['selection'].replace('ligand', self.spec['ligand'])
        self.spec['selected_atoms'] = universe.select_atoms(self.spec['selection'])
        self.spec['selected_atoms'].set_bfactors(1)
        self.spec['refpdb'] = 'omega.pdb'
        self.spec['selected_atoms'].write(os.path.join('input', self.spec['refpdb']))

    def write(self):
        return self.template.render(self.spec)

class ColvarPin:
    _template = """
colvar {
    name {{name}}

    distance {
        group1 {
            dummyatom ( {{ "%8.3f"|format(refx) }}, {{ "%8.3f"|format(refy) }}, {{ "%8.3f"|format(refz) }} )
        }
        group2 {
            atomnumbers {
                {%- for column in selected_atoms|batch(10) %}
                {% for atom in column -%}
                {{ "%6i" | format(atom.number+1) }} 
                {%- endfor -%}
                {% endfor %}
            }
        }
    }
}

harmonic {
    name          {{ name }}
    colvars       {{ name }}
    forceconstant {{ forceconstant|default(100) }}
    centers       {{ centers|default(0.0) }}
}
"""

    def __init__(self, name, universe, **kwargs):
        self.template = Template(self._template)
        self.spec = kwargs
        self.spec['name'] = name
        if self.spec.has_key('receptor'): self.spec['selection'] = self.spec['selection'].replace('receptor', self.spec['receptor'])
        if self.spec.has_key('ligand'): self.spec['selection'] = self.spec['selection'].replace('ligand', self.spec['ligand'])
        self.spec['selected_atoms'] = universe.select_atoms(self.spec['selection'])
        com = self.spec['selected_atoms'].center_of_mass()
        self.spec['refx'], self.spec['refy'], self.spec['refz'] = com

    def write(self):
        return self.template.render(self.spec)



class Colvars(object):
    _template = """colvarstrajfrequency 20
%s
"""
    def __init__(self, conf, psffile, pdbfile):
        self.colvars = []
        self.colvar_keys = {}
        self.conf = conf
        self.universe = mda.Universe(psffile, pdbfile)

    def __iter__(self):
        for colvar in self.colvars:
            yield colvar

    def append(self, colvar_type, name, **kwargs):
        colvar = None
        if colvar_type == 'RMSD':
            colvar = ColvarRMSD(name, self.universe, **kwargs)
        if colvar_type == 'Angle':
            colvar = ColvarAngle(name, self.universe, **kwargs)
        if colvar_type == 'Distance':
            colvar = ColvarDistance(name, self.universe, **kwargs)
        if colvar_type == 'Pin':
            colvar = ColvarPin(name, self.universe, **kwargs)
        if colvar_type == 'Omega':
            colvar = ColvarOmega(name, self.universe, **kwargs)
        if name in self.colvar_keys:
            idx = self.colvar_keys[name]
            self.colvars[idx] = colvar
        else:
            self.colvar_keys[name] = len(self.colvars)
            self.colvars.append(colvar)
        return self

    def write(self, filename):
        colvar_str = ""
        for colvar in self.colvars:
            if colvar.spec['name'] in ['Pin', 'Omega']: continue
            colvar_str += colvar.write() + "\n"

        colvar_str += self.colvars[self.colvar_keys['Pin']].write() + "\n"
        colvar_str += self.colvars[self.colvar_keys['Omega']].write() + "\n"

        fp = open(filename, 'w')
        fp.write(self._template % colvar_str)
        fp.close()
