#!/usr/bin/env python3

import argparse
import numpy
import pandas as pd
import sys

from Bio import SeqIO
from Bio.Seq import Seq


def get_args():
    """
    Define program parameters

    Return
    ------
    argparse object
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
            "genbank",
            metavar="reference-genbank",
            type=str,
            )

    parser.add_argument(
            "-v",
            "--vcf",
            type=str,
            help="Path to vcf file"
            )

    args = parser.parse_args()
    return args


def extract_definitions(json_file):
    """
    CURRENTLY NOT USING
    Extract defining mutations from json file

    Parameters
    ----------
    jons_file : str
        json file containing defining mutation in Constellation format
        See https://github.com/cov-lineages/constellations

    Returns
    -------
    DataFrame of defining mutations
    """
    definitions = pd.read_json(json_file)
    keys = ["lineage", "mutation", "gene", "type", "codon_num", "ref_aa", "alt_aa"]
    var_df = pd.DataFrame(
            [i for j in [[[mutation[key] for key in keys] for mutation in lineage] for lineage in definitions["mutations"]] for i in j],
            columns=keys
            )
    return var_df


def get_cds_from_genbank(genbank_SeqRecord):
    """
    Parameters
    ----------
    genbank_SeqRecord : Bio.SeqRecord.SeqRecord
        Read from GenBank file

    Returns
    -------
    List of cdses

    Note
    ----
    0-indexed
    """
    cdses = [i for i in genbank_SeqRecord.features if i.type == "CDS"]
    gene_names = [i.qualifiers["gene"][0] for i in cdses]
    if len(gene_names) != len(set(gene_names)):
        dup_genes = [i for i in set(gene_names) if gene_names.count(i) > 1]
        gene_str = "gene is" if len(dup_genes) <= 1 else "genes are"
        print("WARNING!! The following "+gene_str+" duplicate in the genbank file. Program may response in weird way.")
        print("\n".join(dup_genes))
    return cdses


def nuc_pos_to_codon_pos(nuc_position, gene_start, codon_start=1, codon_len=3):
    """
    Calculate position of codon on a gene

    Parameters
    ----------
    nuc_position : int
        Position of nucleotide
    gene_start : int
        Position of nucleotide at gene position +1
    codon_start : int
        Number of nucleotide to offset gene position +1
        Default : 1
    codon_len : int
        Length of nucleotide in one codon
        Default : 3

    Returns
    -------
    int
        Position of codon

    Note
    ----
    0-indexed
    """
    return (nuc_position-gene_start+codon_start)//codon_len


def get_codon_base_pos(nuc_position, gene_start, codon_start=1, codon_len=3):
    """
    Calculate base position in a codon

    Parameters
    ----------
    nuc_position : int
        Position of nucleotide
    gene_start : int
        Position of nucleotide at gene position +1
    codon_start : int
        Number of nucleotide to offset gene position +1
        Default : 1
    codon_len : int
        Length of nucleotide in one codon
        Default : 3

    Returns
    -------
    int
        Position of base in a codon

    Note
    ----
    0-indexed
    """
    return (nuc_position-gene_start+codon_start) % codon_len


def codon_pos_to_nuc_pos(codon_pos, gene_start, codon_start=0):
    """
    Calculate nucleotide position

    Parameters
    ----------

    0-indexed
    """
    nuc_pos_start = codon_pos*3
    nuc_pos_start += gene_start
    nuc_pos_start += codon_start
    nuc_pos_end = nuc_pos_start+3
    return nuc_pos_start, nuc_pos_end


def add_alt_codon(alt_codon_dict, position, alt_nuc, frequency, cds):
    """
    Find and add alternative codons to alt_codon_dict

    Parameters
    ----------
    alt_codon_dict : dict
        Nested dict
        {gene: {codon_num : [alt_codon]}}
    position : int
        Position of nucleotide
    alt_nuc : str
        Alternative nucleotide
    frequency : float
        Frequency of alternative nucleotide
    cds : genbank_SeqRecord.features
        with type == "CDS"

    Returns
    -------
    Nested dict with a new alt_codon
    {gene: {codon_num : [alt_codon]}}
    """
    gene = cds.qualifiers["gene"][0]
    codon_num = nuc_pos_to_codon_pos(
            position,
            cds.location.start,
            int(cds.qualifiers["codon_start"][0])-1,
            )
    codon_base_pos = get_codon_base_pos(
            position,
            cds.location.start,
            int(cds.qualifiers["codon_start"][0])-1,
            )
    if gene not in alt_codon_dict.keys():
        alt_codon_dict[gene] = {}
    _pos, codon, _freq = alt_codon_dict[gene].get(codon_num, [0, "xxx", 0.0])
    if codon[codon_base_pos] == "x":
        codon = codon[:codon_base_pos] + alt_nuc + codon[codon_base_pos+1:]
    else:
        # More than one variant exists
        print("Multiple variants not supported")
        print(f"{gene}\t{position}\t{alt_nuc}")
        return alt_codon_dict
    alt_codon_dict[gene][codon_num] = dict()
    alt_codon_dict[gene][codon_num]["nuc"] = position
    alt_codon_dict[gene][codon_num]["codon"] = codon
    alt_codon_dict[gene][codon_num]["frequency"] = frequency

    return alt_codon_dict


def alter_pos(position, cds):
    """
    Calculate positions of nucleotide for cdses with multiple parts

    Parameters
    ----------
    position : int
        Position of nucleotide
    cds : genbank_SeqRecord.features
        with type == "CDS"

    Returns
    -------
    list of int

    Note
    ----
    Too complicated?
    """
    if len(cds.location.parts) > 1:
        positions = []
        if cds.location_operator == "join":
            parts = [part for part in cds.location.parts if part.start <= position < part.end]
            for part in parts:
                alt_position = position
                if part != cds.location.parts[0]:
                    index = cds.location.parts.index(part)
                    for i in range(1, index+1):
                        shift = cds.location.parts[i-1].end - cds.location.parts[i].start
                        alt_position += shift
                positions.append(alt_position)
        else:
            print(f"WARNING!! The operation {cds.location_operator} in this CDS is not support")
        return positions
    return [position]


def make_alt_codon_dict(variant_list, cdses):
    """
    Make a dict containing gene with alternative codons

    parameters
    ----------
    variant_list : list
        List of variants
        [{"CHROM": "Chr1", "POS": 1, "REF": "A", "ID": ".", "ALT", "G", "QUAL": 400, "INFO":{}}, ...]
    cdses : list
        list of genbank_SeqRecord.features with type == "CDS"

    Return
    ------
    Nested dict
    {gene: {codon_num : [alt_codon]}}
    """
    pos_alt_nuc_list = [(var["POS"], var["ALT"], var["INFO"]["AF"]) for var in variant_list]
    alt_codon_dict = {}
    for position, alt_nuc, frequency in pos_alt_nuc_list:
        if len(alt_nuc) != 1:
            continue # Skip indels
        hit_cdses = [cds for cds in cdses if cds.location.start <= position < cds.location.end]
        if not hit_cdses:
            if "nuc" not in alt_codon_dict.keys():
                alt_codon_dict["nuc"] = {}
            alt_codon_dict["nuc"][position] = alt_nuc
            alt_codon_dict["nuc"]["frequency"] = frequency
        else:
            for cds in hit_cdses:
                for alt_position in alter_pos(position, cds):
                    alt_codon_dict = add_alt_codon(alt_codon_dict, alt_position, alt_nuc, float(frequency), cds)
    return alt_codon_dict


def get_nuc_to_codon_coord(cds):
    """
    Calculate nucleotide position of CDS and (codon position and nucleotide position in a codon) coordinates

    Parameters
    ----------
    cds : genbank_SeqRecord.features with type == "CDS"

    Returns
    -------
    dict
    {int: (int, int)}
    {nuc_pos: (codon_pos, nuc_pos_of_codon)}
    {0: (0, 0), 1: (0, 1), 2: (0, 2), 3: (1, 0)}
    """
    nuc_pos = [i for part in cds.location.parts for i in list(range(part.start, part.end))]
    nuc_to_codon = {
            pos: (index//3, index % 3) for index, pos in enumerate([
                    pos for part in cds.location.parts for pos in range(part.start, part.end)
                    ])
            }
    return nuc_to_codon


def get_codon_to_nuc_coord(cds):
    """
    Calculate codon position and nucleotide positions of that codon coordinates


    Parameters
    ----------
    cds : genbank_SeqRecord.features with type == "CDS"

    Returns
    -------
    dict
    {int: (int, int, int)}
    {0: [0. 1, 2], 1: [3, 4, 5]}
    {0: [100, 101, 101], 1: [102, 103, 104]}
    """
    nuc_pos = [i for part in cds.location.parts for i in list(range(part.start, part.end))]
    codon_to_nuc = {
            codon_num: (nuc_pos[codon_num*3], nuc_pos[codon_num*3+1], nuc_pos[codon_num*3+2])
            for codon_num in range(int(len(nuc_pos)/3))
            }
    return codon_to_nuc


def get_alt_codon_seq(alt_codon_dict, cdses, genbank, table=1):
    """

    """
    print_flag = True
    codon_pos_cdses = {cds.qualifiers["gene"][0]: get_codon_to_nuc_coord(cds) for cds in cdses}
    for cds in cdses:
        gene = cds.qualifiers["gene"][0]
        if gene in alt_codon_dict.keys():
            for codon_pos, value in alt_codon_dict[gene].items():
                nuc_pos, alt_codon, freq = value.values()
                ref_codon = "".join(genbank[i] for i in codon_pos_cdses[gene][codon_pos])
                ref_aa = cds.qualifiers["translation"][0]+"*"
                ref_aa = ref_aa[codon_pos]
                alt_codon = "".join([
                        alt_codon[i] if alt_codon[i] != "x"
                        else ref_codon[i]
                        for i in range(3)
                        ])
                alt_aa = str(Seq(alt_codon).translate(table=table))
                alt_aa = alt_aa if alt_aa != ref_aa else "."
                if print_flag:
                    print("\t".join(["pos", "gene", "codon", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "freq"]))
                    print_flag = False
                print(f"{nuc_pos+1}\t{gene}\t{codon_pos+1}\t{ref_codon}\t{alt_codon}\t{ref_aa}\t{alt_aa}\t{freq}")


def classify_mutation(row):
    """
    Give mutation type

    Parameters
    ----------
    row : row from pandas.DataFrame

    Returns
    -------
    str
    """
    len_ref = len(row["ref"])
    len_alt = len(row["alt"])
    if len_ref == 1 and len_alt == 1:
        return "substitution"
    elif len_ref < len_alt:
        return "insertion"
    elif len_ref > len_alt:
        return "deletion"


def import_variant_table(variant_table_file_path):
    """
    Import variant table and convert it to DataFrame

    Parameters
    ----------
    variant_table_file_path : str
        Path to variant table file

    Returns
    -------
    DataFrame

    Notes
    -----
    Convert position to 0-indexed
    """
    variant_df = pd.read_csv(variant_table_file_path, sep="\t", header=None, names=["pos", "ref", "alt"])
    variant_df["pos"] = variant_df.apply(lambda row: int(row["pos"])-1, axis=1)
    variant_df["mutation"] = variant_df.apply(lambda row: classify_mutation(row), axis=1)
    return variant_df


def get_mutation_df(variant_list, cdses):
    """
    Find mutations from alternative nucleotides

    Parameters
    ----------
    variant_list : list
        Alternative nucleotide details
    cdses : list
        List of cdses

    Returns
    -------
    DataFrame
    """
    genes = []
    nuc_pos = []
    mutation_type = []
    codon_pos = []
    pos_in_codon = []
    ref = []
    alt = []

    nuc_to_codon_pos = {
            cds.qualifiers["gene"][0]: get_nuc_to_codon_coord(cds)
            for cds in cdses
            }

    for var in variant_list:
        nuc_number = var["POS"]
        mut_type = "insertion" if len(var["REF"]) < len(var["ALT"]) else "deletion" if len(var["REF"]) > len(var["ALT"]) else "substitution"
        hit_cdses = [cds for cds in cdses if cds.location.start <= nuc_number < cds.location.end]
        if not hit_cdses:
            genes.append("nuc")
            nuc_pos.append(nuc_number)
            mutation_type.append(mut_type)
            codon_pos.append(numpy.nan)
            pos_in_codon.append(numpy.nan)
            ref.append(var["REF"])
            alt.append(var["ALT"])

        else:
            for cds in hit_cdses:
                gene = cds.qualifiers["gene"][0]
                genes.append(gene)
                nuc_pos.append(nuc_number)
                mutation_type.append(mut_type)
                codon_pos.append(nuc_to_codon_pos[gene][nuc_number][0])
                pos_in_codon.append(nuc_to_codon_pos[gene][nuc_number][1])
                ref.append(var["REF"])
                alt.append(var["ALT"])

    mutation_df = pd.DataFrame({
        "nuc_pos": nuc_pos,
        "gene": genes,
        "mutation_type": mutation_type,
        "codon_pos": codon_pos,
        "pos_in_codon": pos_in_codon,
        "ref": ref,
        "alt": alt,
    })

    return mutation_df


def get_mutated_codons(mutation_df):
    """
    Make a nested dict of alternative codons from DataFrame of alternative nucleotide

    Parameters
    ----------
    mutation_df : DataFrame

    Returns
    -------
    dict
    {str: {str: {int: str}}}
    {gene_name: {type_of_variant: [alternative_codon]}}
    """
    mutated_codon_dict = {}
    for row in mutation_df.loc[mutation_df.gene != "nuc"].itertuples():
        nuc_pos = row[1]
        gene = row[2]
        mutation_type = row[3]
        codon_pos = int(row[4])
        pos_in_codon = int(row[5])
        ref = row[6]
        alt = row[7]
        if gene not in mutated_codon_dict.keys():
            mutated_codon_dict[gene] = {}
        if mutation_type not in mutated_codon_dict[gene].keys():
            mutated_codon_dict[gene][mutation_type] = {}

        if mutation_type == "substitution":
            codons = mutated_codon_dict[gene]["substitution"].get(codon_pos, ["xxx"])
            if codons[0][pos_in_codon] != "x":
                new_codons = [codon[:pos_in_codon] + alt + codon[pos_in_codon+1:] for codon in codons]
                codons += new_codons
            else:
                codons = [codon[:pos_in_codon] + alt + codon[pos_in_codon+1:] for codon in codons]
            mutated_codon_dict[gene]["substitution"][codon_pos] = codons

        elif mutation_type == "deletion":
            num_deleted_codons = (len(ref)-1)//3
            frameshift = (len(ref)-1) % 3 != 0
            codon_pos += 1
            mutated_codon_dict[gene]["deletion"][codon_pos] = [nuc_pos, num_deleted_codons, frameshift]

        elif mutation_type == "insertion":
            pass

    return mutated_codon_dict


def format_info(entry):
    """
    Convert INFO column into dict

    Parameters
    ----------
    entry : dict

    Returns
    -------
    dict
    """
    info = entry["INFO"]
    info = info.split(";")
    info = {i: j for i, j in [k.split("=") if "=" in k else [k, k] for k in info]}
    entry["INFO"] = info
    return entry


def parse_vcf(vcf):
    """
    Import vcf file and convert into list of dict

    Parameters
    ----------
    vcf : _io.TextIOWrapper
        vcf file handle

    Returns
    -------
    list

    Notes
    -----
    POS is converted to 0-indexed
    """
    header_section = ""
    variants = []
    for line in vcf:
        if line.startswith("##"):
            header_section += line
        elif line.startswith("#"):
            column_headers = line.strip().replace("#", "")
            column_headers = column_headers.strip().split("\t")
        else:
            line = line.strip().split("\t")
            line = {i: j for i, j in zip(column_headers, line)}
            line["POS"] = int(line["POS"]) - 1
            line = format_info(line)
            variants.append(line)
    return variants


def print_nuc_variant(variants_list):
    """
    Print variant details

    Parameters
    ----------
    variants_list : list
        list of dict from parse_vcf()

    Notes
    -----
    Convert POS back to 1-indexed
    """
    for entry in variants_list:
        print("\t".join([entry["POS"] + 1, entry["REF"], entry["ALT"], f'{float(entry["INFO"]["AF"])*100:>6.2f}']))


def main():
    """
    Main function
    """
    args = get_args()

    # Read vcf file
    if args.vcf:
        with open(args.vcf) as fh:
            variants = parse_vcf(fh)
    else:
        variants = parse_vcf(sys.stdin)

    # Read GenBank file
    genbank = SeqIO.read(args.genbank, "gb")
    cdses = get_cds_from_genbank(genbank)

    # Check amino acid change
    alt_codon_dict = make_alt_codon_dict(variants, cdses)
    get_alt_codon_seq(alt_codon_dict, cdses, genbank, table=1)

    mutation_df = get_mutation_df(variants, cdses)
    mutated_codon_dict = get_mutated_codons(mutation_df)
    deletion_flag = True
    for gene, mutation_types in mutated_codon_dict.items():
        for mutation_type, codons in mutation_types.items():
            if mutation_type == "deletion":
                if deletion_flag:
                    print("DELETION")
                    print("\t".join(["pos", "gene", "codon", "num_codons", "frame_shift"]))
                    deletion_flag = False
                for codon_pos, details in codons.items():
                    print(f"{details[0] + 1}\t{gene}\t{codon_pos+1}\t{details[1]}\t{details[2]}")


if __name__ == "__main__":
    main()
