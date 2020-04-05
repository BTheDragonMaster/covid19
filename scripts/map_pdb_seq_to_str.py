import typing
from pathlib import Path
import prody as pd

def get_sequences_from_fasta_yield(fasta_file: typing.Union[str, Path]) -> tuple:
    """
    Returns (accession, sequence) iterator
    Parameters
    ----------
    fasta_file

    Returns
    -------
    (accession, sequence)
    """
    with open(fasta_file) as f:
        current_sequence = ""
        current_key = None
        for line in f:
            if not len(line.strip()):
                continue
            if "==" in line:
                continue
            if ">" in line:
                if current_key is None:
                    current_key = line.split(">")[1].strip()
                else:
                    yield (current_key, current_sequence)
                    current_sequence = ""
                    current_key = line.split(">")[1].strip()
            else:
                current_sequence += line.strip()
        yield (current_key, current_sequence)


def get_sequences_from_fasta(fasta_file: typing.Union[str, Path]) -> dict:
    """
    Returns dict of accession to sequence from fasta file
    Parameters
    ----------
    fasta_file

    Returns
    -------
    {accession:sequence}
    """
    return {key: sequence for (key, sequence) in get_sequences_from_fasta_yield(fasta_file, prune_headers)}


def get_alpha_indices(protein):
    """
    Get indices of alpha carbons of pd AtomGroup object
    """
    return [i for i, a in enumerate(protein.iterAtoms()) if a.getName() == 'CA']


def map_pdb_sequences_to_structures(pdb_sequences, pdb_folder, map_folder):
    for name_chain in pdb_sequences:
        name, chain = name_chain.split(":")
        pdb_file = pdb_folder / f"{name}.pdb.gz"
        if not pdb_file.exists():
            pdb_file = pdb_folder / f"{name}.cif.gz"
            protein = pd.parseCIF(str(pdb_file), chain=chain.upper())
        else:
            protein = pd.parsePDB(str(pdb_file), chain=chain.upper())
        alpha_indices = get_alpha_indices(protein)
        structure_sequence = protein[alpha_indices].getSequence()
        sequence_file = Path(map_folder) / f"{name_chain}.fasta"
        with open(sequence_file, "w") as f:
            f.write(f">{name}\n{structure_sequence}\n")
            f.write(f">{name}_full\n{pdb_sequences[name_chain]}\n")
        aln_file = Path(map_folder) / f"{name_chain}_aln.fasta"
        
        # REPLACE THIS WITH MUSCLE
        helper.clustal_msa_from_sequences(sequence_file, aln_file, n_iter=0, threads=1)
        
        aln_sequences = helper.get_sequences_from_fasta(aln_file, prune_headers=False)
        pdb_sequence = aln_sequences[name]
        full_sequence = aln_sequences[f"{name}_full"]
        index_1 = 0
        index_2 = 0
        with open(Path(map_folder) / f"{name_chain}.txt", "w") as f:
            f.write("Structure Index\tPDB Index\tSequence Index\n") 
            for i in range(len(pdb_sequence)):
                if pdb_sequence[i] != "-" and full_sequence[i] != "-":
                    f.write(f"{index_1}\t{protein[alpha_indices[index_1]].getResnum()}\t{index_2}\n")
                if pdb_sequence[i] != "-":
                    index_1 += 1
                if full_sequence[i] != "-":
                    index_2 += 1
                    
if __name__ == "__main__":
    sequences = get_sequences_from_fasta("data/pdb_sequences_full.fasta")
    map_pdb_sequences_to_structures(sequences, "data/pdb_files", "data/pdb_residue_mapping")