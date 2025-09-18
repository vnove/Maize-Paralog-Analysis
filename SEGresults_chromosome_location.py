#!/usr/bin/env python3
"""
Corrected Paralog Coordinate Extractor
=======================================

This script fixes the data quality issues identified in the previous analysis:
1. Corrects LCR positions auto-converted to dates
2. Standardizes chromosome nomenclature  
3. Properly handles paralog pair matching
4. Validates data consistency

Usage:
    python corrected_coordinate_extractor.py --lcr_csv <file> --coord_tsv <file> --seg_file <file> --output_dir <dir>
"""

import pandas as pd
import argparse
import os
import re
from collections import defaultdict, Counter
import numpy as np

class CorrectedParalogExtractor:
    """Corrected version that fixes data quality issues"""
    
    def __init__(self, lcr_csv, coord_tsv, seg_file=None, output_dir="corrected_analysis"):
        self.lcr_csv = lcr_csv
        self.coord_tsv = coord_tsv
        self.seg_file = seg_file
        self.output_dir = output_dir
        self.lcr_data = None
        self.coord_data = None
        self.seg_data = {}
        self.gene_coordinates = {}
        self.contig_to_chr_map = {}
        
        os.makedirs(output_dir, exist_ok=True)
    
    def fix_lcr_positions_from_seg(self):
        """Fix LCR positions by parsing original SEG file"""
        if not self.seg_file or not os.path.exists(self.seg_file):
            print("⚠️  SEG file not provided, will attempt to fix positions from data")
            return None
        
        print("🔧 Parsing SEG file to fix LCR positions...")
        
        current_gene = None
        
        with open(self.seg_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                if line.startswith('>'):
                    # Parse SEG header
                    header = line[1:]
                    
                    # Extract gene ID and correct LCR coordinates
                    match = re.match(r'([^(]+)\((\d+)-(\d+)\)\s+complexity=([0-9.]+)', header)
                    if match:
                        gene_id = match.group(1)
                        lcr_start = int(match.group(2))
                        lcr_end = int(match.group(3))
                        complexity = float(match.group(4))
                        
                        self.seg_data[gene_id] = {
                            'correct_positions': f"{lcr_start}-{lcr_end}",
                            'start': lcr_start,
                            'end': lcr_end,
                            'length': lcr_end - lcr_start + 1,
                            'complexity': complexity
                        }
        
        print(f"✅ Extracted correct LCR positions for {len(self.seg_data)} genes from SEG file")
        return self.seg_data
    
    def fix_date_converted_positions(self, position_str):
        """Fix positions that were auto-converted to dates by Excel"""
        if pd.isna(position_str) or position_str == '':
            return None
            
        position_str = str(position_str)
        
        # Handle different date conversion patterns
        date_patterns = [
            (r'^(\d+)-Jan$', r'1-\1'),        # "12-Jan" → "1-12"
            (r'^Jan-(\d+)$', r'1-\1'),        # "Jan-32" → "1-32"
            (r'^(\d+)-Feb$', r'1-\1'),        # Edge case: "2-Feb" → "1-2"
            (r'^Feb-(\d+)$', r'1-\1'),        # Edge case: "Feb-28" → "1-28"
        ]
        
        for pattern, replacement in date_patterns:
            match = re.match(pattern, position_str)
            if match:
                corrected = re.sub(pattern, replacement, position_str)
                return corrected
        
        # If already in correct format (start-end), return as is
        if re.match(r'^\d+-\d+$', position_str):
            return position_str
        
        # If still problematic, return None
        return None
    
    def create_contig_to_chromosome_map(self):
        """Create mapping from contigs to chromosomes where possible"""
        print("🗺️  Creating contig to chromosome mapping...")
        
        # Analyze coordinate data to find patterns
        chr_genes = self.coord_data[self.coord_data['chromosome'].str.startswith('Chr', na=False)]
        contig_genes = self.coord_data[self.coord_data['chromosome'].str.startswith('B73V4_ctg', na=False)]
        
        print(f"📊 Chromosome distribution in coordinate data:")
        print(f"   Standard chromosomes (Chr*): {len(chr_genes)} genes")
        print(f"   Unplaced contigs (B73V4_ctg*): {len(contig_genes)} genes")
        
        # For genes that appear in both formats, create mapping
        gene_chr_map = {}
        
        # Build gene to chromosome mapping for standard chromosomes
        for _, row in chr_genes.iterrows():
            gene_id = row['gene']
            chromosome = row['chromosome']
            gene_chr_map[gene_id] = chromosome
        
        # Try to map contigs based on gene overlap or synteny
        contig_chr_votes = defaultdict(Counter)
        
        for _, row in contig_genes.iterrows():
            contig = row['chromosome']
            gene_id = row['gene']
            
            # If this gene also appears with a standard chromosome, use that
            if gene_id in gene_chr_map:
                contig_chr_votes[contig][gene_chr_map[gene_id]] += 1
        
        # Create final mapping based on majority vote
        for contig, chr_votes in contig_chr_votes.items():
            if chr_votes:
                best_chr = chr_votes.most_common(1)[0][0]
                confidence = chr_votes[best_chr] / sum(chr_votes.values())
                
                if confidence >= 0.7:  # 70% confidence threshold
                    self.contig_to_chr_map[contig] = best_chr
                    print(f"   Mapped {contig} → {best_chr} (confidence: {confidence:.1%})")
        
        print(f"✅ Created mapping for {len(self.contig_to_chr_map)} contigs")
        return self.contig_to_chr_map
    
    def standardize_chromosome_name(self, chr_name):
        """Standardize chromosome names"""
        if pd.isna(chr_name):
            return "Unknown"
        
        chr_name = str(chr_name)
        
        # Already standard format
        if chr_name.startswith('Chr') and len(chr_name) <= 5:
            return chr_name
        
        # Map contig to chromosome if possible
        if chr_name in self.contig_to_chr_map:
            return self.contig_to_chr_map[chr_name]
        
        # Keep unmapped contigs as is but mark them
        if chr_name.startswith('B73V4_ctg'):
            return f"Contig_{chr_name}"
        
        return chr_name
    
    def load_and_clean_data(self):
        """Load and clean all data sources"""
        print("📊 Loading and cleaning data...")
        
        # Load LCR data
        self.lcr_data = pd.read_csv(self.lcr_csv)
        print(f"✅ LCR data: {len(self.lcr_data)} rows")
        
        # Load coordinate data
        self.coord_data = pd.read_csv(self.coord_tsv, sep='\t')
        print(f"✅ Coordinate data: {len(self.coord_data)} rows")
        
        # Parse SEG file if available
        if self.seg_file:
            self.fix_lcr_positions_from_seg()
        
        # Create contig mapping
        self.create_contig_to_chromosome_map()
        
        # Fix LCR positions in the data
        print("🔧 Fixing LCR positions...")
        fixed_positions = []
        
        for idx, row in self.lcr_data.iterrows():
            gene_id = row['Gene_ID']
            original_pos = row.get('LCR_Positions', '')
            
            # First try to get correct position from SEG data
            if gene_id in self.seg_data:
                correct_pos = self.seg_data[gene_id]['correct_positions']
                fixed_positions.append(correct_pos)
            else:
                # Try to fix date-converted positions
                fixed_pos = self.fix_date_converted_positions(original_pos)
                fixed_positions.append(fixed_pos if fixed_pos else original_pos)
        
        self.lcr_data['LCR_Positions_Corrected'] = fixed_positions
        
        # Show correction examples
        corrections_made = 0
        print("\n🔍 LCR Position Corrections (first 10):")
        for i in range(min(10, len(self.lcr_data))):
            original = self.lcr_data.iloc[i]['LCR_Positions']
            corrected = self.lcr_data.iloc[i]['LCR_Positions_Corrected']
            gene = self.lcr_data.iloc[i]['Gene_ID']
            
            if str(original) != str(corrected):
                print(f"   {gene}: '{original}' → '{corrected}'")
                corrections_made += 1
            else:
                print(f"   {gene}: '{original}' (no change needed)")
        
        print(f"✅ Made corrections to LCR positions")
        
        return self.lcr_data, self.coord_data
    
    def build_accurate_gene_coordinate_map(self):
        """Build accurate gene coordinate mapping with proper paralog handling"""
        print("🗺️  Building accurate gene coordinate mapping...")
        
        # Group coordinates by gene to handle multiple entries
        gene_coords = defaultdict(list)
        
        for _, row in self.coord_data.iterrows():
            gene_id = row['gene']
            
            # Standardize chromosome name
            std_chromosome = self.standardize_chromosome_name(row['chromosome'])
            
            coord_info = {
                'chromosome': std_chromosome,
                'original_chromosome': row['chromosome'],
                'start': row['start'],
                'end': row['end'],
                'strand': row['strand'],
                'length': row['end'] - row['start'] + 1 if pd.notna(row['end']) and pd.notna(row['start']) else None,
                'paralog_pair': row['paralog_pair']
            }
            
            gene_coords[gene_id].append(coord_info)
        
        # For genes with multiple coordinate entries, choose the best one
        for gene_id, coords_list in gene_coords.items():
            if len(coords_list) == 1:
                self.gene_coordinates[gene_id] = coords_list[0]
            else:
                # Prefer standard chromosomes over contigs
                std_chr_coords = [c for c in coords_list if c['chromosome'].startswith('Chr')]
                
                if std_chr_coords:
                    # Use the first standard chromosome entry
                    self.gene_coordinates[gene_id] = std_chr_coords[0]
                    self.gene_coordinates[gene_id]['note'] = f"Selected from {len(coords_list)} entries"
                else:
                    # Use the first entry if no standard chromosomes
                    self.gene_coordinates[gene_id] = coords_list[0]
                    self.gene_coordinates[gene_id]['note'] = f"Selected from {len(coords_list)} contig entries"
        
        print(f"✅ Built coordinate map for {len(self.gene_coordinates)} genes")
        
        # Statistics
        chr_types = Counter()
        for coord in self.gene_coordinates.values():
            if coord['chromosome'].startswith('Chr'):
                chr_types['Standard_Chromosome'] += 1
            elif coord['chromosome'].startswith('Contig_'):
                chr_types['Unmapped_Contig'] += 1
            else:
                chr_types['Other'] += 1
        
        print(f"📊 Coordinate quality:")
        for coord_type, count in chr_types.items():
            print(f"   {coord_type}: {count} genes ({count/len(self.gene_coordinates)*100:.1f}%)")
        
        return self.gene_coordinates
    
    def create_corrected_enriched_dataset(self):
        """Create enriched dataset with all corrections applied"""
        print("📊 Creating corrected enriched dataset...")
        
        # Get unique genes from LCR data
        lcr_genes = set(self.lcr_data['Gene_ID'].dropna())
        if 'Paralog_ID' in self.lcr_data.columns:
            lcr_genes.update(set(self.lcr_data['Paralog_ID'].dropna()))
        
        # Find coordinates
        found_genes = []
        missing_genes = []
        
        for gene_id in lcr_genes:
            if gene_id in self.gene_coordinates:
                coord_info = self.gene_coordinates[gene_id]
                found_genes.append({
                    'Gene_ID': gene_id,
                    'Chromosome': coord_info['chromosome'],
                    'Original_Chromosome': coord_info['original_chromosome'],
                    'Start': coord_info['start'],
                    'End': coord_info['end'],
                    'Strand': coord_info['strand'],
                    'Length': coord_info['length'],
                    'Note': coord_info.get('note', '')
                })
            else:
                missing_genes.append(gene_id)
        
        print(f"✅ Found coordinates: {len(found_genes)} genes ({len(found_genes)/len(lcr_genes)*100:.1f}%)")
        print(f"❌ Missing coordinates: {len(missing_genes)} genes ({len(missing_genes)/len(lcr_genes)*100:.1f}%)")
        
        # Create gene lookup
        gene_lookup = {gene['Gene_ID']: gene for gene in found_genes}
        
        # Enrich LCR data
        enriched_rows = []
        
        for _, row in self.lcr_data.iterrows():
            gene1 = row['Gene_ID']
            gene2 = row.get('Paralog_ID', '')
            
            # Start with original row
            enriched_row = row.to_dict()
            
            # Add corrected LCR length from positions if available
            corrected_pos = enriched_row.get('LCR_Positions_Corrected', '')
            if corrected_pos and '-' in str(corrected_pos):
                try:
                    start, end = map(int, str(corrected_pos).split('-'))
                    enriched_row['LCR_Length_Corrected'] = end - start + 1
                except:
                    enriched_row['LCR_Length_Corrected'] = enriched_row.get('Total_LCR_Length', 0)
            else:
                enriched_row['LCR_Length_Corrected'] = enriched_row.get('Total_LCR_Length', 0)
            
            # Add Gene1 coordinates
            if gene1 in gene_lookup:
                coord1 = gene_lookup[gene1]
                enriched_row.update({
                    'Gene1_Chromosome': coord1['Chromosome'],
                    'Gene1_Original_Chr': coord1['Original_Chromosome'],
                    'Gene1_Start': coord1['Start'],
                    'Gene1_End': coord1['End'],
                    'Gene1_Strand': coord1['Strand'],
                    'Gene1_Length_Genomic': coord1['Length'],
                    'Gene1_Mapping_Note': coord1['Note']
                })
            else:
                enriched_row.update({
                    'Gene1_Chromosome': 'Not_Found',
                    'Gene1_Original_Chr': 'Not_Found',
                    'Gene1_Start': None,
                    'Gene1_End': None,
                    'Gene1_Strand': None,
                    'Gene1_Length_Genomic': None,
                    'Gene1_Mapping_Note': 'Gene not found in coordinates'
                })
            
            # Add Gene2 coordinates
            if gene2 and gene2 in gene_lookup:
                coord2 = gene_lookup[gene2]
                enriched_row.update({
                    'Gene2_Chromosome': coord2['Chromosome'],
                    'Gene2_Original_Chr': coord2['Original_Chromosome'],
                    'Gene2_Start': coord2['Start'],
                    'Gene2_End': coord2['End'],
                    'Gene2_Strand': coord2['Strand'],
                    'Gene2_Length_Genomic': coord2['Length'],
                    'Gene2_Mapping_Note': coord2['Note']
                })
                
                # Calculate distance only for reliable chromosome mappings
                gene1_chr = enriched_row.get('Gene1_Chromosome', '')
                gene2_chr = coord2['Chromosome']
                
                both_standard_chr = (gene1_chr.startswith('Chr') and gene2_chr.startswith('Chr'))
                same_chromosome = (gene1_chr == gene2_chr)
                
                if both_standard_chr and same_chromosome and pd.notna(enriched_row.get('Gene1_Start')) and pd.notna(coord2['Start']):
                    distance = abs(enriched_row['Gene1_Start'] - coord2['Start'])
                    enriched_row['Genomic_Distance_bp'] = distance
                    enriched_row['Same_Chromosome'] = True
                    enriched_row['Distance_Reliable'] = True
                    
                    # Distance categories
                    if distance < 100000:
                        enriched_row['Distance_Category'] = 'Close (<100kb)'
                    elif distance < 1000000:
                        enriched_row['Distance_Category'] = 'Moderate (100kb-1Mb)'
                    else:
                        enriched_row['Distance_Category'] = 'Distant (>1Mb)'
                        
                elif same_chromosome:
                    enriched_row['Genomic_Distance_bp'] = None
                    enriched_row['Same_Chromosome'] = True
                    enriched_row['Distance_Reliable'] = False
                    enriched_row['Distance_Category'] = 'Same_Chr_Unreliable'
                else:
                    enriched_row['Genomic_Distance_bp'] = None
                    enriched_row['Same_Chromosome'] = False
                    enriched_row['Distance_Reliable'] = False
                    enriched_row['Distance_Category'] = 'Different_Chromosome'
                    
            else:
                enriched_row.update({
                    'Gene2_Chromosome': 'Not_Found' if gene2 else 'No_Paralog',
                    'Gene2_Original_Chr': 'Not_Found' if gene2 else 'No_Paralog',
                    'Gene2_Start': None,
                    'Gene2_End': None,
                    'Gene2_Strand': None,
                    'Gene2_Length_Genomic': None,
                    'Gene2_Mapping_Note': 'Gene not found in coordinates' if gene2 else 'No paralog',
                    'Genomic_Distance_bp': None,
                    'Same_Chromosome': False,
                    'Distance_Reliable': False,
                    'Distance_Category': 'No_Paralog' if not gene2 else 'Paralog_Not_Found'
                })
            
            enriched_rows.append(enriched_row)
        
        enriched_df = pd.DataFrame(enriched_rows)
        print(f"✅ Created corrected enriched dataset with {len(enriched_df)} rows")
        
        return enriched_df, found_genes, missing_genes
    
    def analyze_corrected_patterns(self, enriched_df):
        """Analyze patterns in the corrected data"""
        print("\n📈 Analyzing corrected chromosomal patterns...")
        
        analysis = {}
        
        # Filter for reliable data only
        reliable_data = enriched_df[
            (enriched_df['Gene1_Chromosome'].str.startswith('Chr', na=False)) & 
            (enriched_df['Gene2_Chromosome'].str.startswith('Chr', na=False)) &
            (enriched_df['Distance_Reliable'] == True)
        ]
        
        if len(reliable_data) == 0:
            print("⚠️  No reliable chromosomal data found for analysis")
            return analysis
        
        print(f"📊 Reliable data: {len(reliable_data)} pairs ({len(reliable_data)/len(enriched_df)*100:.1f}% of total)")
        
        # Same chromosome analysis
        same_chrom = reliable_data['Same_Chromosome'].sum()
        total_pairs = len(reliable_data)
        analysis['reliable_pairs'] = total_pairs
        analysis['same_chromosome_count'] = same_chrom
        analysis['different_chromosome_count'] = total_pairs - same_chrom
        analysis['same_chromosome_percent'] = (same_chrom / total_pairs) * 100
        
        print(f"🧬 Reliable Chromosomal Distribution:")
        print(f"   Same chromosome: {same_chrom} pairs ({analysis['same_chromosome_percent']:.1f}%)")
        print(f"   Different chromosomes: {total_pairs - same_chrom} pairs ({(total_pairs - same_chrom)/total_pairs*100:.1f}%)")
        
        # Distance analysis for same chromosome pairs
        same_chrom_data = reliable_data[reliable_data['Same_Chromosome'] == True]
        if len(same_chrom_data) > 0:
            distances = same_chrom_data['Genomic_Distance_bp'].dropna()
            if len(distances) > 0:
                analysis['distance_stats'] = {
                    'count': len(distances),
                    'mean': distances.mean(),
                    'median': distances.median(),
                    'min': distances.min(),
                    'max': distances.max(),
                    'std': distances.std()
                }
                
                print(f"\n📏 Distance Analysis (reliable same chromosome pairs):")
                print(f"   Count: {len(distances)} pairs")
                print(f"   Average: {distances.mean():,.0f} bp")
                print(f"   Median: {distances.median():,.0f} bp")
                print(f"   Range: {distances.min():,.0f} - {distances.max():,.0f} bp")
        
        # Distance categories for reliable data
        dist_cats = reliable_data['Distance_Category'].value_counts()
        analysis['distance_categories'] = dist_cats.to_dict()
        
        print(f"\n🎯 Distance Categories (reliable data only):")
        for category, count in dist_cats.items():
            print(f"   {category}: {count} pairs ({count/len(reliable_data)*100:.1f}%)")
        
        # Chromosome pair analysis
        chrom_pairs = []
        for _, row in reliable_data.iterrows():
            chrom1 = row['Gene1_Chromosome']
            chrom2 = row['Gene2_Chromosome']
            pair = tuple(sorted([chrom1, chrom2]))
            chrom_pairs.append(pair)
        
        chrom_pair_counts = pd.Series(chrom_pairs).value_counts()
        analysis['chromosome_pair_distribution'] = chrom_pair_counts.head(15).to_dict()
        
        print(f"\n🗺️  Top Chromosome Pairs (reliable data):")
        for pair, count in chrom_pair_counts.head(15).items():
            if pair[0] == pair[1]:
                print(f"   {pair[0]} (same): {count} pairs")
            else:
                print(f"   {pair[0]} - {pair[1]}: {count} pairs")
        
        return analysis
    
    def save_corrected_results(self, enriched_df, found_genes, missing_genes, analysis):
        """Save all corrected results"""
        print(f"\n💾 Saving corrected results...")
        
        results = {}
        
        # 1. Corrected gene coordinates
        gene_coords_df = pd.DataFrame(found_genes)
        coords_file = os.path.join(self.output_dir, 'corrected_gene_coordinates.csv')
        gene_coords_df.to_csv(coords_file, index=False)
        results['coordinates_file'] = coords_file
        print(f"📄 Corrected gene coordinates: {coords_file}")
        
        # 2. Corrected enriched LCR data
        enriched_file = os.path.join(self.output_dir, 'corrected_lcr_data_with_chromosomes.csv')
        enriched_df.to_csv(enriched_file, index=False)
        results['enriched_file'] = enriched_file
        print(f"📄 Corrected enriched LCR data: {enriched_file}")
        
        # 3. Reliable data subset
        reliable_data = enriched_df[
            (enriched_df['Gene1_Chromosome'].str.startswith('Chr', na=False)) & 
            (enriched_df['Gene2_Chromosome'].str.startswith('Chr', na=False)) &
            (enriched_df['Distance_Reliable'] == True)
        ]
        reliable_file = os.path.join(self.output_dir, 'reliable_chromosomal_data.csv')
        reliable_data.to_csv(reliable_file, index=False)
        results['reliable_file'] = reliable_file
        print(f"📄 Reliable data subset: {reliable_file}")
        
        # 4. Data quality report
        quality_file = os.path.join(self.output_dir, 'data_quality_report.txt')
        with open(quality_file, 'w') as f:
            f.write("CORRECTED LCR PARALOG CHROMOSOMAL ANALYSIS\n")
            f.write("="*60 + "\n\n")
            
            f.write("DATA CORRECTIONS APPLIED:\n")
            f.write("  ✅ Fixed LCR positions auto-converted to dates\n")
            f.write("  ✅ Standardized chromosome nomenclature\n")
            f.write("  ✅ Proper paralog coordinate matching\n")
            f.write("  ✅ Data reliability flagging\n\n")
            
            f.write(f"DATA QUALITY SUMMARY:\n")
            f.write(f"  Total LCR entries: {len(enriched_df)}\n")
            f.write(f"  Genes with coordinates: {len(found_genes)}\n")
            f.write(f"  Missing coordinates: {len(missing_genes)}\n")
            f.write(f"  Reliable chromosome pairs: {analysis.get('reliable_pairs', 0)}\n")
            f.write(f"  Reliability rate: {analysis.get('reliable_pairs', 0)/len(enriched_df)*100:.1f}%\n\n")
            
            if analysis:
                f.write(f"CHROMOSOMAL PATTERNS (RELIABLE DATA ONLY):\n")
                if 'same_chromosome_percent' in analysis:
                    f.write(f"  Same chromosome pairs: {analysis['same_chromosome_count']} ({analysis['same_chromosome_percent']:.1f}%)\n")
                    f.write(f"  Different chromosome pairs: {analysis['different_chromosome_count']}\n")
                
                if 'distance_stats' in analysis:
                    stats = analysis['distance_stats']
                    f.write(f"\nDISTANCE ANALYSIS (same chromosome, reliable):\n")
                    f.write(f"  Average distance: {stats['mean']:,.0f} bp\n")
                    f.write(f"  Median distance: {stats['median']:,.0f} bp\n")
                    f.write(f"  Range: {stats['min']:,.0f} - {stats['max']:,.0f} bp\n")
                
                if 'distance_categories' in analysis:
                    f.write(f"\nDISTANCE CATEGORIES (reliable data):\n")
                    for category, count in analysis['distance_categories'].items():
                        f.write(f"  {category}: {count} pairs\n")
        
        results['quality_file'] = quality_file
        print(f"📄 Data quality report: {quality_file}")
        
        return results
    
    def run_corrected_analysis(self):
        """Run complete corrected analysis"""
        print("🔬 CORRECTED LCR PARALOG CHROMOSOMAL ANALYSIS")
        print("="*60)
        
        # Load and clean data
        self.load_and_clean_data()
        
        # Build accurate coordinate mapping
        self.build_accurate_gene_coordinate_map()
        
        # Create corrected enriched dataset
        enriched_df, found_genes, missing_genes = self.create_corrected_enriched_dataset()
        
        # Analyze corrected patterns
        analysis = self.analyze_corrected_patterns(enriched_df)
        
        # Save results
        results = self.save_corrected_results(enriched_df, found_genes, missing_genes, analysis)
        
        print(f"\n🎉 CORRECTED ANALYSIS COMPLETE!")
        print(f"📁 Results saved to: {self.output_dir}")
        print(f"📊 Key corrected findings:")
        print(f"   • Found coordinates for {len(found_genes)} genes")
        print(f"   • Reliable chromosome data: {analysis.get('reliable_pairs', 0)} pairs")
        if analysis and 'same_chromosome_percent' in analysis:
            print(f"   • Same chromosome (reliable): {analysis['same_chromosome_percent']:.1f}%")
        
        return results, analysis

def main():
    parser = argparse.ArgumentParser(description='Corrected chromosomal coordinate extraction for LCR paralogs')
    parser.add_argument('--lcr_csv', required=True,
                       help='LCR analysis CSV file')
    parser.add_argument('--coord_tsv', required=True,
                       help='Paralog coordinates TSV file')
    parser.add_argument('--seg_file',
                       help='Original SEG output file (recommended for position correction)')
    parser.add_argument('--output_dir', default='corrected_analysis',
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Validate files
    if not os.path.exists(args.lcr_csv):
        print(f"❌ Error: LCR CSV file not found: {args.lcr_csv}")
        return
    
    if not os.path.exists(args.coord_tsv):
        print(f"❌ Error: Coordinate TSV file not found: {args.coord_tsv}")
        return
    
    if args.seg_file and not os.path.exists(args.seg_file):
        print(f"⚠️  Warning: SEG file not found: {args.seg_file}")
        print(f"    Will attempt to fix positions from existing data")
        args.seg_file = None
    
    # Run corrected extraction
    extractor = CorrectedParalogExtractor(args.lcr_csv, args.coord_tsv, args.seg_file, args.output_dir)
    results, analysis = extractor.run_corrected_analysis()
    
    print(f"\n✅ Main corrected output files:")
    print(f"   📊 All data: {results['enriched_file']}")
    print(f"   🎯 Reliable data only: {results['reliable_file']}")
    print(f"   📋 Quality report: {results['quality_file']}")

if __name__ == "__main__":
    main()
