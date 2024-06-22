"""Module for visualizing protein-protein interactions using PyMOL"""
import itertools
import os
from pathlib import Path
from typing import Iterable

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
        """Python PyMOL wrapper to visualize protein-protein interactions. The wrapper is based on 
        a more general wrapper implemented in the `Graphein package <https://graphein.ai/>`_.

        Args:
            port (int, optional): Port to use for communication with PyMOL. Defaults to PyMOLPORT
                from Graphein.
        """
        self.port = port
        self.session_buffer = []

    def display_ppi(
        self,
        ppi_path: Path,
        reuse_session: bool = False,
        colors: Iterable[str] = ('hotpink', 'greencyan'),
        residue_color_sets: Iterable[str] = ('reds+magentas', 'greens+yellows+oranges'),
        swap_colors: bool = False,
        transparency: float = 0.95,
        color_by_residues: bool = False,
        sticks: bool = False,
        letters: bool = True
    ):
        """Display protein-protein interaction using PyMOL. If used in a Jupyter notebook, the 
        method may also be used to show a static image of the interaction.

        Args:
            ppi_path (Path): Path to a .pdb file containing a PPI. It is recommended to use
                a file produced by the ``ppiref.extraction.PPIExtractor`` class.
            reuse_session (bool, optional): If set to True, displays PPI in the same session as 
                during the previous call. Otherwise, creates a new PyMOL sessions. Defaults to False.
            colors (Iterable[str], optional): Two colors to use for two interacting proteins. 
                Defaults to ``('hotpink', 'greencyan')``.
            residue_color_sets (Iterable[str], optional): If ``color_by_residue`` is set to True,
                the two provided palettes will be used to color each type of residue in a different
                color in two interacting proteins. Please see ``PYMOL_COLOR_SETS`` from the same
                modile for the list of available options. Defaults to ``('reds+magentas',
                'greens+yellows+oranges')``.
            swap_colors (bool, optional): Swap colors for two interacting proteins. Defaults to False.
            transparency (float, optional): Transparency factor from the [0, 1] range. Defaults to 0.95.
            color_by_residues (bool, optional): Color residues of different amino acid types with
                different colors. Defaults to False.
            sticks (bool, optional): Show amino acids in the stick representation (more detailed).
                Defaults to False.
            letters (bool, optional): Show one-letter codes of amino acid types. Defaults to True.
        """
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
    
        # Create new session if necessary
        while pymol is None:
            try:
                pymol = MolViewer(port=self.port)
            except:
                pass
            self.port += 1

        # Swap colors for protein partners
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
        if os.path.exists(f'{pdb_id}.cif'):
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
