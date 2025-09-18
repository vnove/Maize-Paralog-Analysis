#!/usr/bin/env python3
"""
SEG Output Parser for Maize Paralog LCR Analysis
Extracts comprehensive LCR data from SEG output file
Usage: python seg_parser.py input_seg_file.out output_lcr_summary.csv
"""

import sys
import csv
import re
from collections import defaultdict

def parse_seg_output(seg_file_path):
    """
    Parse SEG output file and extract LCR information for each gene
    Returns: Dictionary with gene-centric LCR data
    """
    gene_data = defaultdict(lambda: {
        'gene_id': '',
        'lcr_regions': [],
        'total_lcr_length': 0,
        'lcr_count': 0,
        'avg_complexity': 0.0,
        'max_complexity': 0.0,
        'min_complexity': float('inf'),
        'lcr_positions': []
    })
    
    print(f"Processing SEG file: {seg_file_path}")
    
    with open(seg_file_path, 'r') as file:
        lines = file.readlines()
    
    i = 0
    total_genes_processed = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Check if line is a header (starts with >)
        if line.startswith('>'):
            # Extract gene ID using regex
            gene_match = re.search(r'>(Zm\d+d\d+_T\d+)', line)
            if gene_match:
                gene_id = gene_match.group(1)
                
                # Extract position information
                pos_match = re.search(r'\((\d+)-(\d+)\)', line)
                if pos_match:
                    start_pos = int(pos_match.group(1))
                    end_pos = int(pos_match.group(2))
                    lcr_length = end_pos - start_pos + 1
                    
                    # Extract complexity score
                    complexity_match = re.search(r'complexity=([\d.]+)', line)
                    complexity = float(complexity_match.group(1)) if complexity_match else 0.0
                    
                    # Get the sequence (next line)
                    sequence = ''
                    if i + 1 < len(lines) and not lines[i + 1].startswith('>'):
                        sequence = lines[i + 1].strip()
                    
                    # Store LCR region data
                    lcr_region = {
                        'start': start_pos,
                        'end': end_pos,
                        'length': lcr_length,
                        'complexity': complexity,
                        'sequence': sequence
                    }
                    
                    # Add to gene data
                    if gene_data[gene_id]['gene_id'] == '':
                        gene_data[gene_id]['gene_id'] = gene_id
                    
                    gene_data[gene_id]['lcr_regions'].append(lcr_region)
                    gene_data[gene_id]['total_lcr_length'] += lcr_length
                    gene_data[gene_id]['lcr_count'] += 1
                    gene_data[gene_id]['lcr_positions'].append(f"{start_pos}-{end_pos}")
                    
                    # Update complexity statistics
                    if complexity > gene_data[gene_id]['max_complexity']:
                        gene_data[gene_id]['max_complexity'] = complexity
                    if complexity < gene_data[gene_id]['min_complexity']:
                        gene_data[gene_id]['min_complexity'] = complexity
        
        i += 1
    
    # Calculate average complexity for each gene
    for gene_id in gene_data:
        if gene_data[gene_id]['lcr_count'] > 0:
            total_complexity = sum(region['complexity'] for region in gene_data[gene_id]['lcr_regions'])
            gene_data[gene_id]['avg_complexity'] = total_complexity / gene_data[gene_id]['lcr_count']
            total_genes_processed += 1
        else:
            gene_data[gene_id]['min_complexity'] = 0.0
    
    print(f"Processed {total_genes_processed} genes with LCRs")
    print(f"Total LCR regions found: {sum(data['lcr_count'] for data in gene_data.values())}")
    
    return dict(gene_data)

def write_lcr_summary(gene_data, output_file):
    """
    Write gene-centric LCR summary to CSV file
    """
    print(f"Writing LCR summary to: {output_file}")
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = [
            'Gene_ID',
            'LCR_Count',
            'Total_LCR_Length',
            'Avg_Complexity',
            'Max_Complexity',
            'Min_Complexity',
            'LCR_Positions',
            'First_LCR_Start',
            'Last_LCR_End',
            'LCR_Span',
            'Has_Multiple_LCRs'
        ]
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        # Sort genes by ID for consistent output
        sorted_genes = sorted(gene_data.keys())
        
        for gene_id in sorted_genes:
            data = gene_data[gene_id]
            
            # Calculate additional metrics
            first_lcr_start = min(region['start'] for region in data['lcr_regions']) if data['lcr_regions'] else 0
            last_lcr_end = max(region['end'] for region in data['lcr_regions']) if data['lcr_regions'] else 0
            lcr_span = last_lcr_end - first_lcr_start + 1 if data['lcr_regions'] else 0
            
            row = {
                'Gene_ID': gene_id,
                'LCR_Count': data['lcr_count'],
                'Total_LCR_Length': data['total_lcr_length'],
                'Avg_Complexity': round(data['avg_complexity'], 3),
                'Max_Complexity': round(data['max_complexity'], 3),
                'Min_Complexity': round(data['min_complexity'], 3) if data['min_complexity'] != float('inf') else 0,
                'LCR_Positions': ';'.join(data['lcr_positions']),
                'First_LCR_Start': first_lcr_start,
                'Last_LCR_End': last_lcr_end,
                'LCR_Span': lcr_span,
                'Has_Multiple_LCRs': 'Yes' if data['lcr_count'] > 1 else 'No'
            }
            
            writer.writerow(row)
    
    print(f"LCR summary written successfully!")

def print_statistics(gene_data):
    """
    Print summary statistics about the LCR data
    """
    print("\n" + "="*50)
    print("LCR ANALYSIS STATISTICS")
    print("="*50)
    
    lcr_counts = [data['lcr_count'] for data in gene_data.values()]
    lcr_count_dist = defaultdict(int)
    for count in lcr_counts:
        lcr_count_dist[count] += 1
    
    print(f"Total genes analyzed: {len(gene_data)}")
    print(f"Total LCR regions: {sum(lcr_counts)}")
    print(f"Average LCRs per gene: {sum(lcr_counts)/len(gene_data):.2f}")
    
    print("\nLCR Count Distribution:")
    for lcr_count in sorted(lcr_count_dist.keys()):
        gene_count = lcr_count_dist[lcr_count]
        percentage = (gene_count / len(gene_data)) * 100
        print(f"  {lcr_count} LCRs: {gene_count} genes ({percentage:.1f}%)")
    
    # Top genes with most LCRs
    top_genes = sorted(gene_data.items(), key=lambda x: x[1]['lcr_count'], reverse=True)[:10]
    print(f"\nTop 10 genes with most LCRs:")
    for gene_id, data in top_genes:
        print(f"  {gene_id}: {data['lcr_count']} LCRs, {data['total_lcr_length']} bp total")
    
    print("="*50)

def main():
    if len(sys.argv) != 3:
        print("Usage: python seg_parser.py <input_seg_file> <output_csv_file>")
        print("Example: python seg_parser.py valid_paralog_nucleotide_sequences_seg_W12_T1_E1.3.out lcr_summary.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        # Parse SEG output
        gene_data = parse_seg_output(input_file)
        
        # Print statistics
        print_statistics(gene_data)
        
        # Write summary CSV
        write_lcr_summary(gene_data, output_file)
        
        print(f"\n✅ Processing complete!")
        print(f"📊 Input: {input_file}")
        print(f"📁 Output: {output_file}")
        print(f"🧬 Ready for LCR-paralogy correlation analysis!")
        
    except FileNotFoundError:
        print(f"❌ Error: Could not find input file '{input_file}'")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
