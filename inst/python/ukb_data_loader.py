#!/usr/bin/env python3
"""
UK Biobank RAP Data Downloader (Spark-based)

Author: Nan He, Southern Medical University, Basic Medical Sciences, Department of Bioinformatics

Usage:
    python ukb_data_loader.py demographic --ids 31,53,21022 --output population.csv
    python ukb_data_loader.py demographic --id-file field_ids.txt --output population.csv
    python ukb_data_loader.py metabolites --output metabolites.csv
    python ukb_data_loader.py metabolites --non-ratio --output metabolites_non_ratio.csv
"""

import pyspark
import dxpy
import dxdata
import pandas as pd
from packaging.version import Version, InvalidVersion
import re
import argparse
import os


# Global variables
spark = None
clinical_participant_dataset = None


def init_spark():
    """Initialize Spark session and load dataset."""
    global spark, clinical_participant_dataset
    
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    
    # Load dataset
    dispensed_dataset_id = dxpy.find_one_data_object(
        typename="Dataset", 
        name="app*.dataset", 
        folder="/", 
        name_mode="glob"
    )
    dataset_id = dispensed_dataset_id['id']
    clinical_dataset = dxdata.load_dataset(id=dataset_id)
    clinical_participant_dataset = clinical_dataset['participant']
    
    print(f"Connected to dataset: {dataset_id}")


def fields_for_id(field_id):
    """Get all fields for a given UKB field ID, sorted by version."""
    field_id = str(field_id)
    fields = clinical_participant_dataset.find_fields(
        name_regex=r'^p{}(_i\d+)?(_a\d+)?$'.format(field_id)
    )
    
    def safe_version(name):
        version_pattern = r'^p{}\D*(\d+(\.\d+)*)?$'.format(field_id)
        match = re.match(version_pattern, name)
        if match and match.group(1):
            try:
                return Version(match.group(1))
            except InvalidVersion:
                pass
        return Version('0.0.0')
    
    return sorted(fields, key=lambda f: safe_version(f.name))


def field_names_for_id(field_id):
    """Get field names for a given UKB field ID."""
    return [f.name for f in fields_for_id(field_id)]


def load_ids_from_file(filepath):
    """
    Load field IDs from a text file.
    Supports: one ID per line, comma-separated, or space-separated.
    Lines starting with # are ignored.
    """
    ids = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Support comma or space separated
            parts = re.split(r'[,\s]+', line)
            ids.extend([p.strip() for p in parts if p.strip()])
    return ids


def download_demographic(field_ids, output_path):
    """
    Download demographic data by field IDs.
    
    Args:
        field_ids: List of UKB field IDs (e.g., ['31', '53', '21022'])
        output_path: Output CSV file path
    """
    print(f"Downloading demographic data for {len(field_ids)} fields...")
    
    # Get all field names with prefix
    field_names = sum([field_names_for_id(fid) for fid in field_ids], [])
    
    # Add eid
    all_fields = ['eid'] + field_names
    print(f"Total columns: {len(all_fields)}")
    
    # Retrieve and save
    df = clinical_participant_dataset.retrieve_fields(
        names=all_fields, 
        engine=dxdata.connect()
    )
    df.toPandas().to_csv(output_path, index=False)
    print(f"Saved to: {output_path}")


def download_metabolites(output_path, non_ratio=False):
    """
    Download NMR metabolites data.
    
    Args:
        output_path: Output CSV file path
        non_ratio: If True, download only non-ratio metabolites from reference file
    """
    if non_ratio:
        print("Downloading non-ratio metabolites data...")
        # Load field IDs from reference file
        script_dir = os.path.dirname(os.path.abspath(__file__))
        ref_file = os.path.join(script_dir, '..', 'extdata', 'metabolites_non_ratio.txt')
        
        # Read tab-separated file
        ref_df = pd.read_csv(ref_file, sep='\t')
        all_fields = ['eid'] + ref_df['meta_ID'].tolist()
    else:
        print("Downloading all metabolites data...")
        # NMR fields: 20280-20281 and 23400-23648
        nmr_fields_1 = [f"p{i}_i0" for i in range(20280, 20282)]
        nmr_fields_2 = [f"p{i}_i0" for i in range(23400, 23649)]
        all_fields = ['eid'] + nmr_fields_1 + nmr_fields_2
    
    print(f"Total columns: {len(all_fields)}")
    
    # Retrieve and save
    df = clinical_participant_dataset.retrieve_fields(
        names=all_fields, 
        engine=dxdata.connect()
    )
    df.toPandas().to_csv(output_path, index=False)
    print(f"Saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='UK Biobank RAP Data Downloader')
    subparsers = parser.add_subparsers(dest='command', help='Data type to download')
    
    # Demographic subcommand
    demo_parser = subparsers.add_parser('demographic', help='Download demographic data')
    demo_parser.add_argument('--ids', type=str, default=None,
                            help='Comma-separated UKB field IDs (e.g., 31,53,21022)')
    demo_parser.add_argument('--id-file', type=str, default=None,
                            help='File containing UKB field IDs (one per line or comma-separated)')
    demo_parser.add_argument('-o', '--output', type=str, default='population.csv',
                            help='Output file path')
    
    # Metabolites subcommand
    meta_parser = subparsers.add_parser('metabolites', help='Download metabolites data')
    meta_parser.add_argument('--non-ratio', action='store_true',
                            help='Download only non-ratio metabolites (168 fields)')
    meta_parser.add_argument('-o', '--output', type=str, default='metabolites.csv',
                            help='Output file path')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Initialize Spark
    init_spark()
    
    if args.command == 'demographic':
        # Load IDs from file or command line
        if args.id_file:
            field_ids = load_ids_from_file(args.id_file)
        elif args.ids:
            field_ids = [x.strip() for x in args.ids.split(',')]
        else:
            print("Error: Either --ids or --id-file is required")
            return 1
        download_demographic(field_ids, args.output)
    elif args.command == 'metabolites':
        download_metabolites(args.output, non_ratio=args.non_ratio)
    
    print("Done!")
    return 0


if __name__ == "__main__":
    exit(main())
