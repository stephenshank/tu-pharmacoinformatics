import argparse
import csv
from sys import argv

from rdkit.Chem import MolToSmiles
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.AllChem import FragmentOnBonds
from rdkit.Chem.rdmolops import GetMolFrags


def clean_smiles(smiles):
    bracket_index = smiles.find(']')
    return smiles[bracket_index + 1: ]


def ionize_a_fragment(mol):
    mol_weight = MolWt(mol)
    mol_smiles = clean_smiles(MolToSmiles(mol))
    hydrogen = MolFromSmiles('[H]')
    hydrogen_weight = MolWt(hydrogen)
    sodium = MolFromSmiles('[Na]')
    sodium_weight = MolWt(sodium)
    fragment_data = [
        {
            'smiles': mol_smiles,
            'descriptor': 'M',
            'weight': mol_weight,
            'mz': mol_weight
        },
        {
            'smiles': mol_smiles,
            'descriptor': 'M+H',
            'weight': mol_weight + hydrogen_weight,
            'mz': mol_weight + hydrogen_weight 
        },
        {
            'smiles': mol_smiles,
            'descriptor': 'M+Na',
            'weight': mol_weight + sodium_weight,
            'mz': mol_weight + sodium_weight
        }
    ]
    return fragment_data


def fragment_and_ionize(smiles, bonds):
    mol = MolFromSmiles(smiles)
    all_ions = ionize_a_fragment(mol)
    for bond in bonds:
        fragments = GetMolFrags(FragmentOnBonds(mol, [bond]), asMols=True)
        for fragment in fragments:
            all_ions.extend(ionize_a_fragment(fragment))
    return all_ions


DESCRIPTION = """
fragment - Fragment and ionize a molecule in-silico
By Stephen D. Shank, Ph. D.
sshank@temple.edu
Acme Computational Molecular Evolution Group
http://lab.hyphy.org/
"""


def fragment_cli():
    "Full command line interface function."
    if len(argv) == 1:
        print(DESCRIPTION)
        print("For more information: regal --help")
        sys.exit(0)
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-s',
        '--smiles',
        help='SMILES representation of molecule',
        dest='smiles',
        required=True
    )
    parser.add_argument(
        '-b',
        '--bonds',
        help='Bonds to break (comma delimited)',
        dest='bonds',
        required=True
    )
    parser.add_argument(
        '-o',
        '--output',
        help='output (CSV)',
        dest='output',
        required=True
    )
    args = parser.parse_args()
    bonds = [int(i) for i in args.bonds.split(',')]
    all_fragments_and_ions = fragment_and_ionize(args.smiles, bonds)
    with open(args.output, 'w') as output_file:
        fieldnames = all_fragments_and_ions[0].keys()
        csv_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        csv_writer.writeheader()
        for fragmented_ion in all_fragments_and_ions:
            csv_writer.writerow(fragmented_ion)


if __name__ == '__main__':
    fragment_cli()
