#!/usr/bin/env python3
"""
Extract cluster information from HDBSCAN clustering results.

This script parses all CL_XX_info.json and CL_XX_isolates.txt files in a given folder
and creates a comprehensive dictionary with cluster information.
"""

import argparse
import json
import yaml
import os
import re
from pathlib import Path
import numpy as np


def parse_cluster_files(cl_folder):
    """
    Parse all cluster files in the given folder.
    
    Args:
        cl_folder (str): Path to the cluster folder
        
    Returns:
        dict: Dictionary with cluster information
    """
    cl_folder = Path(cl_folder)
    clusters = {}
    
    # Find all CL_XX_info.json files
    info_files = list(cl_folder.glob("CL_*_info.json"))
    
    for info_file in info_files:
        # Extract cluster number from filename
        match = re.search(r'CL_(\d+)_info\.json', info_file.name)
        if not match:
            continue
            
        cluster_num = match.group(1)
        cluster_name = f"CL_{cluster_num}"
        
        # Read info file
        with open(info_file, 'r') as f:
            info_data = json.load(f)
        
        # Read corresponding isolates file
        isolates_file = cl_folder / f"CL_{cluster_num}_isolates.txt"
        if isolates_file.exists():
            isolates = np.loadtxt(isolates_file, dtype=str).tolist()
            # Handle case where there's only one isolate (numpy returns a scalar)
            if isinstance(isolates, str):
                isolates = [isolates]
        else:
            isolates = []
        
        if len(isolates) == 0:
            raise ValueError(f"No isolates found for cluster {cluster_name}. Check if the isolates file {isolates_file} exists and is not empty.")

        # Create cluster entry
        clusters[cluster_name] = {
            'isolates': isolates,
            'info': info_data,
            'size': len(isolates)
        }
    
    return clusters


def get_species_name(cl_folder):
    """
    Extract species name from the cluster folder path.
    
    Args:
        cl_folder (str): Path to the cluster folder
        
    Returns:
        str: Species name
    """
    cl_folder = Path(cl_folder)
    # The species name should be the parent directory of the cluster type folder
    # e.g., clusters/kpneumoniae/hdbscan_clusters -> kpneumoniae
    return cl_folder.parent.name


def parse_args():
    parser = argparse.ArgumentParser(description='Extract cluster information from HDBSCAN clustering results')
    parser.add_argument('--cl-folder', required=True, 
                        help='Path to the cluster folder (e.g., clusters/kpneumoniae/hdbscan_clusters)')
    parser.add_argument('--cl-json', required=True,
                        help='Output JSON file path')
    parser.add_argument('--cl-yaml', required=True,
                        help='Output YAML file path')
    
    return parser.parse_args()


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()

    # Parse cluster files
    print(f"Parsing cluster files in: {args.cl_folder}")
    clusters = parse_cluster_files(args.cl_folder)
    
    # Get species name
    species_name = get_species_name(args.cl_folder)
    
    # Save JSON output
    print(f"Saving JSON output to: {args.cl_json}")
    with open(args.cl_json, 'w') as f:
        json.dump(clusters, f, indent=2)
    
    # Prepare YAML output
    cluster_names = sorted(clusters.keys(), key=lambda name: clusters[name]['size'], reverse=True)
    yaml_data = {species_name: cluster_names}
    
    # Save YAML output
    print(f"Saving YAML output to: {args.cl_yaml}")
    with open(args.cl_yaml, 'w') as f:
        yaml.dump(yaml_data, f, default_flow_style=False)
    
    print(f"Extracted {len(clusters)} clusters for species '{species_name}'")
    print(f"Total isolates across all clusters: {sum(cluster['size'] for cluster in clusters.values())}")