import copy
import itertools
import os
from pathlib import Path

import Bio.PDB
from graphein.utils.pymol import PORT as PyMOLPORT
from graphein.utils.pymol import MolViewer

from ppiref.utils.residue import BASE_AMINO_ACIDS_3
from ppiref.utils.ppipath import path_to_pdb_id


PYMOL_COLOR_SETS = {
   'reds': [
        'red',
        'tv_red',
        'raspberry',
        'darksalmon',
        'salmon',
        'deepsalmon',
        'warmpink',
        'firebrick',
        'ruby',
        'chocolate',
        'brown',
    ],
    'magentas': [
        'magenta',
        'lightmagenta',
        'hotpink',
        'pink',
        'lightpink',
        'dirtyviolet',
        'violet',
        'violetpurple',
        'purple',
        'deeppurple'
    ],
    'greens': [
        'green',
        'tv_green',
        'chartreuse',
        'splitpea',
        'smudge',
        'palegreen',
        'limegreen',
        'lime',
        'limon',
        'forest'
    ],
    'yellows': [
        'yellow',
        'tv_yellow',
        'paleyellow',
        'yelloworange',
        'limon',
        'wheat',
        'sand'
    ],
    'oranges': [
        'orange',
        'tv_orange',
        'brightorange',
        'lightorange',
        'yelloworange',
        'olive',
        'deepolive'
    ]
}


class PyMOL:
    def __init__(
        self,
        port: int = PyMOLPORT
    ):
        self.port = port
        self.session_buffer = []

    def display_ppi(
        self,
        ppi_path: Path,
        reuse_session=False,
        colors=('hotpink', 'greencyan'),
        residue_color_sets=('reds+magentas', 'greens+yellows+oranges'),
        swap_colors=False,
        transparency=0.95,
        color_by_residues=False,
        sticks=False,
        letters=True
    ):
        # Init PyMOL session
        if type(reuse_session) == int:
            pymol = self.session_buffer[reuse_session]
        elif reuse_session is True:
            if len(self.session_buffer) > 0:
                pymol = self.session_buffer[-1]
            else:
                pymol = None
        elif reuse_session is False:
            pymol = None
        else:
            raise ValueError('Wrong `reuse_session` value.')
    
        # Create new if necessary
        while pymol is None:
            try:
                pymol = MolViewer(port=self.port)
            except:
                pass
            self.port += 1

        if swap_colors:
            colors = list(reversed(colors))
            residue_color_sets = list(reversed(residue_color_sets))

        # Prepare session for visualiztion
        pymol.start()
        pymol.delete('all')
        pdb_id = path_to_pdb_id(ppi_path)
        pymol.fetch(f'{pdb_id}')

        # Remove other
        pymol.do('remove solvent')

        # Construct necessary selections
        model = Bio.PDB.PDBParser().get_structure('', ppi_path)[0]
        chain_sels = []
        iface_sels = []
        for chain in model:
            # Add chain selection
            chain_sels.append(f'(chain {chain.get_id()})')

            # Add chain interface selection
            sel = [f'resi {res.get_id()[1]}' for res in chain]
            sel = ' or '.join(sel)
            sel = f'(chain {chain.get_id()} and ({sel}))'
            iface_sels.append(sel)

        # Set transparency style
        pymol.set('transparency_mode', 3)
        pymol.set('sphere_transparency', transparency, 'all')
        pymol.set('cartoon_transparency', transparency, 'all')
        pymol.set('stick_transparency', transparency, 'all')

        # Set cartoon style
        pymol.do('set cartoon_flat_sheets, 0')
        pymol.do('set cartoon_smooth_loops, 0')
        pymol.do('set ray_trace_mode, 3')

        # Set colors
        for sel in iface_sels:
            pymol.set('cartoon_transparency', 0.0, sel)
            pymol.set('stick_transparency', 0.0, sel)
        pymol.color(colors[0], chain_sels[0])
        pymol.color(colors[1], chain_sels[1])

        if color_by_residues:
            self._color_by_residues(pymol, residue_color_sets[0], iface_sels[0])
            self._color_by_residues(pymol, residue_color_sets[1], iface_sels[1])

        # Add residue labels
        if letters:
            for sel in iface_sels:
                pymol.do(f'label n. ca and {sel}, one_letter[resn]')
                pymol.do(f'set label_color, grey, {sel}')
            # pymol.do('set label_size, -0.75')
            pymol.do('set label_size, -1.50')
            pymol.do('set label_font_id, 10')

        # Show sticks on the interface
        if sticks:
            for sel in iface_sels:
                pymol.do(f'show sticks, {sel}')
            pymol.do('color atomic, (not elem C)')

        # Zoom in
        pymol.do(f'zoom ({iface_sels[0]} or {iface_sels[1]})')

        # Remove fog
        pymol.do('set fog, off')

        # Store PyMOL object to keep window open
        self.session_buffer.append(pymol)

        # Delete fetched PDB file
        os.remove(f'{pdb_id}.cif')

        return pymol.display()
    
    def _color_by_residues(self, pymol, colorset_name, selection='all'):
        for aa in BASE_AMINO_ACIDS_3:
            pymol.select(aa,f'resn {aa}')
        pymol.select('none')

        colors = self._get_color_set(colorset_name)
        code = dict(zip(BASE_AMINO_ACIDS_3, colors))

        for aa in code:
            line='color ' + code[aa]+',' + aa + '&' + selection
            pymol.do(line)

        # pymol.do('set ray_trace_mode, 1')

    def _get_color_set(self, name):
        colors = [PYMOL_COLOR_SETS[n] for n in name.split('+')]
        return list(itertools.chain.from_iterable(colors))
