# MAVISp - various utilities for MAVISp
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

three_to_one = { 'ALA' : 'A',
                 'CYS' : 'C',
                 'ASP' : 'D',
                 'GLU' : 'E',
                 'PHE' : 'F',
                 'GLY' : 'G',
                 'HIS' : 'H',
                 'ILE' : 'I',
                 'LYS' : 'K',
                 'LEU' : 'L',
                 'MET' : 'M',
                 'ASN' : 'N',
                 'PRO' : 'P',
                 'GLN' : 'Q',
                 'ARG' : 'R',
                 'SER' : 'S',
                 'THR' : 'T',
                 'VAL' : 'V',
                 'TRP' : 'W',
                 'TYR' : 'Y'}

one_to_three_hgvsp = {'A': 'Ala',
                      'C': 'Cys',
                      'D': 'Asp',
                      'E': 'Glu',
                      'F': 'Phe',
                      'G': 'Gly',
                      'H': 'His',
                      'I': 'Ile',
                      'K': 'Lys',
                      'L': 'Leu',
                      'M': 'Met',
                      'N': 'Asn',
                      'P': 'Pro',
                      'Q': 'Gln',
                      'R': 'Arg',
                      'S': 'Ser',
                      'T': 'Thr',
                      'V': 'Val',
                      'W': 'Trp',
                      'Y': 'Tyr'}

three_to_one_hgvsp = { v: k for k, v in one_to_three_hgvsp.items() }

def mutation_to_HGVSp(mutation):
    return f"p.{one_to_three_hgvsp[mutation[0]]}{mutation[1:-1]}{one_to_three_hgvsp[mutation[-1]]}"

