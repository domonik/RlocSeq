import sqlite3

import pandas as pd
import gffutils
from tempfile import TemporaryDirectory
import os

GFF_COLNAMES = ['seqid', 'source', 'featuretype', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

FIVE_PRIME_NAMES = ["5'UTR", "5UTR", "five_prime_UTR"]
THREE_PRIME_NAMES = ["3'UTR", "3UTR", "three_prime_UTR"]

def parse_attributes(attributes):
    """
    Parse the attributes field of a GFF3 file and return a dictionary of key-value pairs.

    :param attributes: The attributes string from a GFF3 file.
    :return: A dictionary of key-value pairs.
    """
    if pd.isna(attributes):
        return {}
    return dict(item.split('=') for item in attributes.split(';') if '=' in item)


def read_gff3(filepath):
    """
    Read a GFF3 file into a pandas DataFrame with appropriate column names.

    :param filepath: Path to the GFF3 file.
    :return: A pandas DataFrame containing the GFF3 data.
    """
    column_names = GFF_COLNAMES
    df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=column_names)
    attributes_df = df['attributes'].apply(parse_attributes).apply(pd.Series)

    # Combine the original DataFrame with the expanded attributes DataFrame
    df = pd.concat([df.drop(columns=['attributes']), attributes_df], axis=1)
    return df


def filter_gff_by_interval(df, seqid, start, end):
    """
    Filter the GFF3 DataFrame to include only rows where the coordinates overlap with the specified interval.

    :param df: The GFF3 DataFrame.
    :param start: The start of the interval.
    :param end: The end of the interval.
    :return: A filtered pandas DataFrame.
    """
    # Filter based on the overlap condition

    condition = (df['start'] <= end) & (df['end'] >= start)
    df = df[df["seqid"] == seqid]
    filtered_df = df[condition]

    return filtered_df


def find_gene_coordinates(gff, gene, anno_type):
    df = gff[gff[anno_type].str.contains(gene, case=False) == True]
    return df


def coords_from_gff_row(row):
    return row["seqid"], row["start"], row.end, row.strand


ID_SPEC = {
    'gene': ["ID", 'gene_id', 'geneID'],
    'transcript': ["ID", 'transcript_id', 'transcriptID'],
}


def read_gff_via_gffutils(file, dbname=None, fix=True):
    if dbname is None:
        dir = TemporaryDirectory()
        os.path.join(dir.name, "GFF.db")
    if os.path.exists(dbname):
        db = gffutils.FeatureDB(dbfn=dbname)
    else:
        db = gffutils.create_db(file, dbfn=dbname, force=True, keep_order=True, merge_strategy="merge", id_spec=ID_SPEC)
        if fix:
            fix_gff_db(dbname)
    return db


def intervals_overlap(i1, i2):
    return i1.start <= i2.end and i2.start <= i1.end


def nearest_distance(i1, i2):
    # Case 1: If i1 is completely before i2
    if i1.end < i2.start:
        return i2.start - i1.end
    # Case 2: If i2 is completely before i1
    elif i2.end < i1.start:
        return i1.start - i2.end
    # Case 3: If the intervals overlap
    else:
        return 0



def optimize_gff_lines(gff: gffutils.FeatureDB):
    genes = gff.execute(
        """
        SELECT f.*
        FROM features f
        LEFT JOIN relations r ON f.id = r.child
        WHERE f.featuretype = 'gene' OR r.parent IS NULL
        ORDER BY f.seqid, f.start, f.strand;
        """
    )
    lines = []
    mapping = {}
    attributes = {}
    old_seqid = None
    for row in genes:
        seq_id = row["seqid"]
        if seq_id != old_seqid:
            lines = []
            old_seqid = seq_id
        idx = row["id"]
        if idx in mapping:
            raise KeyError("Key in mapping already exists")
        start, end = row['start'], row['end']
        placed = False
        keys = eval(row["attributes"]).keys()
        for key in keys:
            if key in attributes:
                attributes[key] += 1
            else: attributes[key] = 0


        # Try to place the span in an existing line
        for i, line_end in enumerate(lines):
            if start > line_end:
                mapping[idx] = i
                lines[i] = end
                placed = True
                break
        # If no suitable line was found, create a new line
        if not placed:
            mapping[idx] = len(lines)
            lines.append(end)
    sorted_keys = sorted(attributes, key=attributes.get, reverse=True)
    return mapping, sorted_keys


def fix_gff_db(dbfn: str):
    db = gffutils.FeatureDB(dbfn)
    n_genes = db.count_features_of_type("gene")
    if n_genes == 0:
        print("No genes detected - creating genes for CDS")

    has_transcripts = bool(db.count_features_of_type("mRNA"))
    last_gene_or_transcript = None
    print("checking entries")
    relations = []
    for entry in db.all_features():
        if entry.featuretype == "gene" or entry.featuretype == "mRNA":
            last_gene_or_transcript = entry

        elif entry.featuretype in [*FIVE_PRIME_NAMES, *THREE_PRIME_NAMES, "CDS"]:
            parents = list(db.parents(entry.id))
            if not parents:
                if intervals_overlap(last_gene_or_transcript, entry):
                    relations.append((last_gene_or_transcript, [entry], False))
                else:
                    print("Non valid gff3 entry outside gene annotation detected. try to fix")
                    res = list(db.region(strand=entry.strand, start=entry.start -10, end=entry.end+10, seqid=entry.seqid, featuretype="gene"))
                    if res:
                        nearest = min(res, key=lambda r: nearest_distance(r, entry))
                        new_children = list(db.children(nearest)) + [entry]

                        if nearest.end < entry.start:
                            print("feature before. Extending gene downstream")
                            if entry.featuretype == "CDS":
                                continue
                            db.delete([nearest])
                            nearest.end = entry.end
                            delete = True

                        elif nearest.start > entry.end:
                            if entry.featuretype == "CDS":
                                continue
                            db.delete([nearest])
                            nearest.start = entry.start
                            delete = True
                            print("feature after extending gene upstream")
                        else:
                            delete = False
                            print("overlapping")

                        relations.append((nearest, new_children, delete))
    del db
    db = gffutils.FeatureDB(dbfn)
    for gene, targets, need_to_add in relations:
        if need_to_add:
            db.update([gene], keep_order=True, merge_strategy="merge", id_spec="ID")
        listed_children = [child.id for child in db.children(gene.id)]
        for target in targets:
            if target.id not in listed_children:
                db.add_relation(gene, target, level=1)
















