#!/usr/bin/env python

# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

import argparse
import os
import subprocess
import sys

import pandas as pd

def load_cufflinks_fpkm(fpkm_files):
    sample_gene_expressions = []
    for (i, f) in enumerate(fpkm_files):
        data = pd.read_csv(f, sep='\t')
        data['SampleId'] = i
        sample_gene_expressions.append(data)
    sample_gene_expression = pd.concat(sample_gene_expressions)

    # Filter to OK FPKM counts
    sample_gene_expression = sample_gene_expression[sample_gene_expression['FPKM_status'] == 'OK']
    return sample_gene_expression


def main(args):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cufflinks-files",
        nargs='+',
        required=True,
        help="List of cufflinks FPKM files to use"
    )

    # parser.add_argument(
    #     "--cufflinks-dir",
    #     help="Directory which contains cufflinks FPKM files"
    # )

    parser.add_argument(
        "--gep-output",
        default=".tmp.gep.out",
        help="Output file for intermediate GEP file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output file final results"
    )

    parser.add_argument(
        "--cibersort-output",
        default=".tmp.cibersort.out",
        help="Output file for intermediate CIBERSORT file"
    )

    parser.add_argument(
        "--mixture-file",
        required=True,
        help="Path to mixture file to use with CIBERSORT (i.e. LM22.txt)"
    )

    parser.add_argument(
        "--cibersort-jar",
        required=True,
        help="Path to CIBERSORT jar (i.e. CIBERSORT.jar)"
    )

    parser.add_argument(
        "--min-fpkm",
        type=int,
        default=5,
        help="Minimum FPKM for gene in samples (speeds up CIBERSORT)"
    )

    parser.add_argument(
        "--print-top-n",
        action='store_true',
        default=False,
        help="Print top N cell types per sample"
    )

    parser.add_argument(
        "--n",
        type=int,
        default=3,
        help="Print top N cell types per sample"
    )

    args = parser.parse_args()

    if args.cufflinks_files:
        fpkm_files = args.cufflinks_files
    elif args.cufflinks_dir:
        fpkm_files = os.list_dir(args.cufflinks_dir)
        fpkm_files = [f for f in fpkm_files if 'fpkm' in f]

    gene_fpkm = load_cufflinks_fpkm(fpkm_files)

    # pivot gene_fpkm
    gene_fpkm_pivot = pd.pivot_table(gene_fpkm,
                      index=['SampleId'],
                      values=['FPKM'],
                      columns=['gene_short_name'],
                      fill_value=0)
    gene_fpkm_pivot.columns = gene_fpkm_pivot.columns.get_level_values(1)

    # Drop infrequently occurring genes
    gene_fpkm_pivot.drop(
        labels=gene_fpkm_pivot.columns[gene_fpkm_pivot.sum() < args.min_fpkm], 
        inplace=True,
        axis=1
    )

    gene_fpkm_pivot.T.to_csv(args.gep_output, sep='\t')

    # Start R Server
    subprocess.check_call(
        ["R", "CMD", "Rserve", "--no-save"],
        #stdout=open("/dev/null"),
        #stderr="/dev/null",
    )

    # Run CIBERSORT
    # java -jar CIBERSORT.jar -M Mixture -B Signature_Matrix [Options]
    subprocess.check_call(
        ["java", 
        "-jar", 
        args.cibersort_jar, 
        "-B", 
        args.mixture_file,
        "-M",
        args.gep_output],
        stdout=open(args.cibersort_output, 'w')
    )

    cibersort_data = pd.read_csv(args.cibersort_output, sep='\t', comment='>')
    cibersort_data.drop(labels=['Column'], inplace=True, axis=1)
    cibersort_data.set_index(pd.Series(data=fpkm_files, name='SampleId'), inplace=True)
    cibersort_data.to_csv(args.output, index=True, sep='\t')

    if args.print_top_n:
        non_cell_type_columns = cibersort_data.columns[-4:]
        cibersort_cell_types = cibersort_data.drop(labels=non_cell_type_columns, axis=1).reset_index()
        cibersort_cell_types_melted = \
            pd.melt(cibersort_cell_types, id_vars=['SampleId'])
        print(cibersort_cell_types_melted.sort(
                    ['SampleId', 'value'], 
                    ascending=[1, 0]
            ).groupby(['SampleId']).head(args.n))


if __name__ == '__main__':
    main(sys.argv)