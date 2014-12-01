from jinja2 import Template
import MDAnalysis as mda

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
        refpositionfile input/{{refpdb}}
        refpositioncol  B
    }
}

harmonic {
    name          {{ name }}
    colvars       {{ name }}
    forceconstant {{ forceconstant|default(10.0) }}
    centers       {{ min|default(0.0) }}
}
"""

    def __init__(self, name, universe, **kwargs):
        self.template = Template(self._template)
        self.spec = kwargs
        self.spec['name'] = name
        self.spec['selected_atoms'] = universe.selectAtoms(self.spec['selection'])
        print kwargs

    def write(self):
        return self.template.render(self.spec)

class ColvarAngle:
    _template = """
colvar {
    name {{name}}

    {{ angletype }}
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
        self.spec['refatoms'] = [universe.selectAtoms(self.spec['refatoms'][ref]) for ref in self.spec['angle']]
        self.spec['angletype'] = 'angle' if len(self.spec['refatoms']) == 3 else 'dihedral'
        print kwargs

    def write(self):
        return self.template.render(self.spec)

class ColvarDistance:
    _template = """
colvar {
    name {{name}}

    distance
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
        self.spec['refatoms'] = [universe.selectAtoms(self.spec['refatoms'][ref]) for ref in self.spec['distance']]
        print kwargs

    def write(self):
        return self.template.render(self.spec)


class Colvars:
    _template = """colvarstrajfrequency 20
%s
"""

    def __init__(self, conf):
        self.colvars = []
        self.conf = conf
        self.universe = mda.Universe(conf['psffile'], conf['pdbfile'])

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
        self.colvars.append(colvar)
        return self

    def write(self, filename):
        colvar_str = ""
        for colvar in self.colvars:
            colvar_str += colvar.write() + "\n"

        fp = open(filename, 'w')
        fp.write(self._template % colvar_str)
        fp.close()
