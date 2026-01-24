#!/usr/bin/env python3
"""
UK Biobank RAP Proteomics Data Downloader

Author: Nan He, Southern Medical University, Basic Medical Sciences, Department of Bioinformatics

Usage:
    python main.py                      # Download all proteins to current directory
    python main.py --output ./data      # Download to specified directory
    python main.py --batch-size 100     # Custom batch size
"""

import pandas as pd
import subprocess
import dxpy
import os
import glob
import time
import argparse


def get_dataset_id():
    """Get UK Biobank dataset ID from DNAnexus platform."""
    dispensed_dataset_id = dxpy.find_one_data_object(
        typename='Dataset', 
        name='app*.dataset', 
        folder='/', 
        name_mode='glob'
    )['id']
    project_id = dxpy.find_one_project()["id"]
    dataset = f"{project_id}:{dispensed_dataset_id}"
    return dataset


def download_data_dictionary(dataset, output_dir):
    """Download data dictionary to output directory."""
    original_dir = os.getcwd()
    os.chdir(output_dir)
    try:
        cmd = ["dx", "extract_dataset", dataset, "-ddd", "--delimiter", ","]
        subprocess.check_call(cmd)
        data_dict_csv = glob.glob(os.path.join(output_dir, "*.data_dictionary.csv"))[0]
        return data_dict_csv
    finally:
        os.chdir(original_dir)


def get_protein_fields(data_dict_path, entity="olink_instance_0"):
    """Extract protein field names from data dictionary."""
    data_dict_df = pd.read_csv(data_dict_path)
    field_names = list(
        data_dict_df.loc[data_dict_df["entity"] == entity, "name"].values
    )
    field_names_str = [f"{entity}.{f}" for f in field_names]
    return field_names_str


def download_proteins(dataset, field_batches, output_dir, batch_size, delay=2.0):
    """Download protein data in batches."""
    print("Starting batch extraction...")
    print(f"Total {len(field_batches)} batches, max {batch_size} fields per batch")
    print("=" * 50)
    
    successful_files = []
    failed_batches = []
    
    for i, batch in enumerate(field_batches):
        try:
            field_names_batch = ",".join(batch)
            output_filename = f"protein_batch_{i+1}.csv"
            output_path = os.path.join(output_dir, output_filename)
            
            print(f"Processing batch {i+1}/{len(field_batches)}...")
            print(f"  Fields: {len(batch)}")
            print(f"  Output: {output_filename}")
            
            cmd = [
                'dx',
                'extract_dataset',
                dataset,
                '--fields',
                field_names_batch,
                '--delimiter',
                ',',
                '--output',
                output_path
            ]
            
            start_time = time.time()
            subprocess.check_call(cmd)
            end_time = time.time()
            
            successful_files.append(output_filename)
            print(f"  Success! Time: {end_time - start_time:.2f}s")
            
            if i < len(field_batches) - 1:
                print(f"  Waiting {delay}s...")
                time.sleep(delay)
                
        except subprocess.CalledProcessError as e:
            failed_batches.append((i + 1, str(e)))
            print(f"  Failed: {e}")
        except Exception as e:
            failed_batches.append((i + 1, str(e)))
            print(f"  Error: {e}")
        
        print("-" * 30)
    
    # Summary
    print("\n" + "=" * 50)
    print("Summary:")
    print(f"Success: {len(successful_files)} files")
    if failed_batches:
        print(f"Failed: {len(failed_batches)} batches")
        for batch_num, error in failed_batches:
            print(f"  - Batch {batch_num}: {error}")
    print("=" * 50)
    
    return successful_files, failed_batches


def merge_files(successful_files, output_dir, id_column="olink_instance_0.eid"):
    """Merge all batch files into one."""
    print("\nMerging batch files...")
    
    all_dataframes = []
    for filename in successful_files:
        file_path = os.path.join(output_dir, filename)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            all_dataframes.append(df)
            print(f"  {filename}: {df.shape[1]} cols, {len(df)} rows")
    
    if not all_dataframes:
        print("No files to merge")
        return None
    
    # Merge by eid column
    if id_column in all_dataframes[0].columns:
        print(f"\nMerging by {id_column}...")
        merged_df = all_dataframes[0]
        for i, df in enumerate(all_dataframes[1:], 1):
            merged_df = merged_df.merge(df, on=id_column, how='outer')
            print(f"  Merged file {i+1}... Shape: {merged_df.shape}")
    else:
        print("\nConcatenating horizontally...")
        merged_df = pd.concat(all_dataframes, axis=1)
    
    # Save merged file
    merged_path = os.path.join(output_dir, "protein_all_merged.csv")
    merged_df.to_csv(merged_path, index=False)
    
    print(f"\nMerge complete!")
    print(f"  Final shape: {merged_df.shape}")
    print(f"  Saved to: protein_all_merged.csv")
    
    return merged_df


def print_file_sizes(file_list, output_dir):
    """Print file size statistics."""
    print("\nFile sizes:")
    for filename in file_list:
        file_path = os.path.join(output_dir, filename)
        if os.path.exists(file_path):
            size_mb = os.path.getsize(file_path) / (1024 * 1024)
            print(f"  {filename}: {size_mb:.2f} MB")


def main():
    parser = argparse.ArgumentParser(description='UK Biobank RAP Proteomics Data Downloader')
    parser.add_argument('-o', '--output', type=str, default='.', help='Output directory')
    parser.add_argument('--batch-size', type=int, default=200, help='Fields per batch')
    parser.add_argument('--delay', type=float, default=2.0, help='Delay between batches (seconds)')
    parser.add_argument('--no-merge', action='store_true', help='Skip merging batch files')
    args = parser.parse_args()
    
    output_dir = args.output
    batch_size = args.batch_size
    os.makedirs(output_dir, exist_ok=True)
    
    print("UK Biobank RAP Proteomics Data Downloader")
    print("Author: Nan He, Southern Medical University")
    print("=" * 50)
    
    # Step 1: Get dataset ID
    print("\nGetting dataset ID...")
    dataset = get_dataset_id()
    print(f"Dataset: {dataset}")
    
    # Step 2: Download data dictionary
    print("\nDownloading data dictionary...")
    data_dict_path = download_data_dictionary(dataset, output_dir)
    print(f"Data dictionary: {data_dict_path}")
    
    # Step 3: Get protein fields
    print("\nExtracting protein fields...")
    field_names = get_protein_fields(data_dict_path)
    print(f"Total fields: {len(field_names)}")
    
    # Step 4: Create batches
    field_batches = []
    for i in range(0, len(field_names), batch_size):
        batch = field_names[i:i + batch_size]
        field_batches.append(batch)
    
    print(f"\nSplit into {len(field_batches)} batches")
    for i, batch in enumerate(field_batches):
        print(f"  Batch {i+1}: {len(batch)} fields")
    
    # Step 5: Download protein data
    print()
    successful_files, failed_batches = download_proteins(
        dataset, field_batches, output_dir, batch_size, args.delay
    )
    
    # Step 6: Merge files
    if not args.no_merge and successful_files:
        merged_df = merge_files(successful_files, output_dir)
        
        # Print file sizes
        all_files = successful_files + ['protein_all_merged.csv']
        print_file_sizes(all_files, output_dir)
    
    print(f"\nDone! Generated {len(successful_files)} batch files")
    print(f"Output directory: {output_dir}")
    
    return 0 if not failed_batches else 1


if __name__ == "__main__":
    exit(main())
