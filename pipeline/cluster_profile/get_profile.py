from argparse import ArgumentParser
from collections import defaultdict
import os, shutil
from seqcluster.libs.read import map_to_precursors


def get_cluster(in_file, keep):
    seqs = {}
    loci = []
    cluster_name = ""
    with open(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("C "):
                if cluster_name == keep:
                    return seqs, loci
                seqs = {}
                loci = []
                cluster_name = line.strip().split(" ")[1]
            if line.startswith("L "):
                loci = [keep] + line.strip().split(" ")[2:]
            if line.startswith("S "):
                s, f = line.strip().split(" ")[1:]
                seqs[s] = f


if __name__ == "__main__":
    parser = ArgumentParser(description="Check installation")
    parser.add_argument("--json", required=True, help="json file.")
    parser.add_argument("--keep", required=True, help="cluster id.")
    parser.add_argument("--ref", help="Path to aligner index.")
    parser.add_argument("--out", help="output file.")

    args = parser.parse_args()
    os.mkdir(args.out)

    seqs, loci = get_cluster(args.json, args.keep)
    map_to_precursors(seqs.keys(), seqs.values(), {args.keep: [loci]}, args)
