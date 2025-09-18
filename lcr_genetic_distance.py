def question2_genetic_distances(self):
        """
        Question 2: Calculate genetic distances of genes with vs without LCRs
        to compare their evolutionary rates
        """
        print("\n=== QUESTION 2: Genetic Distance Analysis ===")
        
        # Use comprehensive gene information
        if self.all_genes_info is None:
            print("No comprehensive gene information available. Using only LCR dataset.")
            return {}
            
        genes_with_lcrs = [g for g, info in self.all_genes_info.items() if info['has_lcr']]
        genes_without_lcrs = [g for g, info in self.all_genes_info.items() if not info['has_lcr']]
        
        print(f"Total unique genes: {len(self.all_genes_info)}")
        print(f"Genes with LCRs: {len(genes_with_lcrs)}")
        print(f"Genes without LCRs: {len(genes_without_lcrs)}")
        
        # Method 1: Use paralog pairs from LCR dataset (for genes with LCRs)
        def calculate_lcr_gene_distances():
            distances = []
            for idx, row in self.df.iterrows():
                if (row['Has_Paralog_Bool'] and 
                    row['Distance_Reliable'] and
                    pd.notna(row['Genomic_Distance_bp']) and
                    row['Gene_ID'] in genes_with_lcrs):
                    
                    normalized_distance = row['Genomic_Distance_bp'] / row['Gene1_Length_Genomic'] if row['Gene1_Length_Genomic'] > 0 else 0
                    distances.append({
                        'gene_id': row['Gene_ID'],
                        'paralog_id': row['Paralog_ID'],
                        'genomic_distance': row['Genomic_Distance_bp'],
                        'normalized_distance': normalized_distance,
                        'same_chromosome': row['Same_Chromosome'],
                        'gene_length': row['Gene1_Length_Genomic'],
                        'lcr_count': row['Number_of_LCRs'],
                        'gene_group': 'with_LCRs'
                    })
            return distances
        
        # Method 2: Find paralog pairs from paralog_coordinates.tsv for genes without LCRs
        def calculate_no_lcr_gene_distances():
            distances = []
            
            if self.paralog_coords_df is not None:
                print("Processing paralog pairs from paralog_coordinates.tsv...")
                print(f"Paralog coordinates file has {len(self.paralog_coords_df)} rows")
                
                # Debug: Check what genes we're looking for
                print(f"Looking for {len(genes_without_lcrs)} genes without LCRs")
                print(f"Sample genes without LCRs: {list(genes_without_lcrs)[:5]}")
                
                # Debug: Check what genes are in paralog coordinates
                unique_genes_in_coords = set(self.paralog_coords_df['gene'].unique())
                print(f"Unique genes in paralog_coordinates.tsv: {len(unique_genes_in_coords)}")
                print(f"Sample genes in coordinates: {list(unique_genes_in_coords)[:5]}")
                
                # Check overlap
                overlap = set(genes_without_lcrs) & unique_genes_in_coords
                print(f"Genes without LCRs that are in paralog coordinates: {len(overlap)}")
                print(f"Sample overlap: {list(overlap)[:5]}")
                
                # Group by paralog_pair to get both genes in each pair
                paralog_groups = self.paralog_coords_df.groupby('paralog_pair')
                print(f"Total paralog pairs in file: {len(paralog_groups)}")
                
                pairs_with_no_lcr_genes = 0
                
                for pair_name, group in paralog_groups:
                    if len(group) == 2:  # Should have exactly 2 genes per pair
                        gene1_row = group.iloc[0]
                        gene2_row = group.iloc[1]
                        
                        gene1_id = gene1_row['gene']
                        gene2_id = gene2_row['gene']
                        
                        # Check if EITHER gene is in our "without LCRs" group (not both required)
                        gene1_no_lcr = gene1_id in genes_without_lcrs
                        gene2_no_lcr = gene2_id in genes_without_lcrs
                        
                        # Also check if either gene has LCRs
                        gene1_has_lcr = gene1_id in genes_with_lcrs
                        gene2_has_lcr = gene2_id in genes_with_lcrs
                        
                        if gene1_no_lcr or gene2_no_lcr:
                            pairs_with_no_lcr_genes += 1
                            
                            # Debug: Print first few pairs found
                            if pairs_with_no_lcr_genes <= 5:
                                print(f"Found pair {pairs_with_no_lcr_genes}: {gene1_id} (LCR: {gene1_has_lcr}) ↔ {gene2_id} (LCR: {gene2_has_lcr})")
                            
                            # Get coordinates
                            try:
                                gene1_chr = gene1_row['chromosome']
                                gene1_start = int(gene1_row['start'])
                                gene1_end = int(gene1_row['end'])
                                gene1_length = abs(gene1_end - gene1_start)
                                
                                gene2_chr = gene2_row['chromosome']
                                gene2_start = int(gene2_row['start'])
                                gene2_end = int(gene2_row['end'])
                                gene2_length = abs(gene2_end - gene2_start)
                                
                                # Calculate genomic distance
                                if gene1_chr == gene2_chr:
                                    # Same chromosome - calculate actual distance
                                    genomic_distance = abs(gene1_start - gene2_start)
                                    same_chromosome = True
                                    
                                    # Normalize by gene length
                                    normalized_distance_1 = genomic_distance / gene1_length if gene1_length > 0 else 0
                                    normalized_distance_2 = genomic_distance / gene2_length if gene2_length > 0 else 0
                                else:
                                    # Different chromosomes - still record but with different handling
                                    genomic_distance = 1000000000  # Large number for inter-chromosomal
                                    same_chromosome = False
                                    normalized_distance_1 = 0  # Can't meaningfully normalize inter-chromosomal
                                    normalized_distance_2 = 0
                                
                                # Add entries for genes without LCRs only
                                if gene1_no_lcr:
                                    distances.append({
                                        'gene_id': gene1_id,
                                        'paralog_id': gene2_id,
                                        'genomic_distance': genomic_distance,
                                        'normalized_distance': normalized_distance_1,
                                        'same_chromosome': same_chromosome,
                                        'gene_length': gene1_length,
                                        'paralog_length': gene2_length,
                                        'gene_chr': gene1_chr,
                                        'paralog_chr': gene2_chr,
                                        'lcr_count': 0,
                                        'gene_group': 'without_LCRs',
                                        'paralog_has_lcr': gene2_has_lcr
                                    })
                                
                                if gene2_no_lcr:
                                    distances.append({
                                        'gene_id': gene2_id,
                                        'paralog_id': gene1_id,
                                        'genomic_distance': genomic_distance,
                                        'normalized_distance': normalized_distance_2,
                                        'same_chromosome': same_chromosome,
                                        'gene_length': gene2_length,
                                        'paralog_length': gene1_length,
                                        'gene_chr': gene2_chr,
                                        'paralog_chr': gene1_chr,
                                        'lcr_count': 0,
                                        'gene_group': 'without_LCRs',
                                        'paralog_has_lcr': gene1_has_lcr
                                    })
                                    
                            except (ValueError, KeyError) as e:
                                print(f"Error processing pair {pair_name}: {e}")
                                continue
                
                print(f"Paralog pairs involving genes without LCRs: {pairs_with_no_lcr_genes}")
                print(f"Total distance measurements for genes without LCRs: {len(distances)}")
                
                if len(distances) == 0:
                    print("WARNING: No distances calculated for genes without LCRs!")
                    print("This suggests gene IDs don't match between FASTA and paralog_coordinates files")
                    
                    # Debug: Show some examples
                    print("\nDebugging gene ID matching:")
                    print("FASTA gene IDs (first 5):", list(self.sequence_genes)[:5])
                    print("Paralog coords gene IDs (first 5):", list(unique_genes_in_coords)[:5])
                    print("LCR dataset gene IDs (first 5):", list(self.df['Gene_ID'])[:5])
                
            return distances
        
        # Calculate distances for both groups
        lcr_distances = calculate_lcr_gene_distances()
        no_lcr_distances = calculate_no_lcr_gene_distances()
        
        print(f"Found {len(lcr_distances)} paralog pairs for genes with LCRs")
        print(f"Found {len(no_lcr_distances)} paralog pairs for genes without LCRs")
        
        # Extract numerical values for statistical analysis
        lcr_distance_values = [d['normalized_distance'] for d in lcr_distances if d['normalized_distance'] > 0]
        no_lcr_distance_values = [d['normalized_distance'] for d in no_lcr_distances if d['normalized_distance'] > 0]
        
        lcr_genomic_values = [d['genomic_distance'] for d in lcr_distances if d['genomic_distance'] > 0]
        no_lcr_genomic_values = [d['genomic_distance'] for d in no_lcr_distances if d['genomic_distance'] > 0]
        
        # Calculate gene length statistics using comprehensive gene info
        lcr_gene_lengths = [info['gene_length'] for gene_id, info in self.all_genes_info.items() 
                          if info['has_lcr'] and info['gene_length'] > 0]
        no_lcr_gene_lengths = [info['gene_length'] for gene_id, info in self.all_genes_info.items() 
                             if not info['has_lcr'] and info['gene_length'] > 0]
        
        # Statistical comparisons
        distance_stats = None
        genomic_distance_stats = None
        length_stats = None
        
        if len(lcr_distance_values) > 0 and len(no_lcr_distance_values) > 0:
            distance_stats = mannwhitneyu(lcr_distance_values, no_lcr_distance_values, alternative='two-sided')
            
        if len(lcr_genomic_values) > 0 and len(no_lcr_genomic_values) > 0:
            genomic_distance_stats = mannwhitneyu(lcr_genomic_values, no_lcr_genomic_values, alternative='two-sided')
        
        if len(lcr_gene_lengths) > 0 and len(no_lcr_gene_lengths) > 0:
            length_stats = mannwhitneyu(lcr_gene_lengths, no_lcr_gene_lengths, alternative='two-sided')
        
        # Correlation analysis: LCR count vs evolutionary distance
        lcr_distance_correlation = None
        if len(lcr_distances) > 1:
            lcr_counts = [d['lcr_count'] for d in lcr_distances]
            distances = [d['normalized_distance'] for d in lcr_distances]
            # Filter out zero distances for correlation
            valid_pairs = [(lc, dist) for lc, dist in zip(lcr_counts, distances) if dist > 0]
            if len(valid_pairs) > 1:
                lcr_counts_valid = [x[0] for x in valid_pairs]
                distances_valid = [x[1] for x in valid_pairs]
                lcr_distance_correlation = pearsonr(lcr_counts_valid, distances_valid)
        
        results = {
            'total_genes_analyzed': len(self.all_genes_info),
            'genes_with_lcrs_count': len(genes_with_lcrs),
            'genes_without_lcrs_count': len(genes_without_lcrs),
            'lcr_paralog_pairs_analyzed': len(lcr_distances),
            'no_lcr_paralog_pairs_analyzed': len(no_lcr_distances),
            'mean_normalized_distance_with_lcrs': np.mean(lcr_distance_values) if lcr_distance_values else 0,
            'mean_normalized_distance_without_lcrs': np.mean(no_lcr_distance_values) if no_lcr_distance_values else 0,
            'median_normalized_distance_with_lcrs': np.median(lcr_distance_values) if lcr_distance_values else 0,
            'median_normalized_distance_without_lcrs': np.median(no_lcr_distance_values) if no_lcr_distance_values else 0,
            'mean_genomic_distance_with_lcrs': np.mean(lcr_genomic_values) if lcr_genomic_values else 0,
            'mean_genomic_distance_without_lcrs': np.mean(no_lcr_genomic_values) if no_lcr_genomic_values else 0,
            'mean_gene_length_with_lcrs': np.mean(lcr_gene_lengths) if lcr_gene_lengths else 0,
            'mean_gene_length_without_lcrs': np.mean(no_lcr_gene_lengths) if no_lcr_gene_lengths else 0,
            'normalized_distance_statistical_test': distance_stats,
            'genomic_distance_statistical_test': genomic_distance_stats,
            'gene_length_statistical_test': length_stats,
            'lcr_count_distance_correlation': lcr_distance_correlation,
            'detailed_lcr_distances': lcr_distances,
            'detailed_no_lcr_distances': no_lcr_distances
        }
        
        self.results['question2'] = results
        
        print(f"Paralog pairs analyzed - with LCRs: {len(lcr_distances)}")
        print(f"Paralog pairs analyzed - without LCRs: {len(no_lcr_distances)}")
        print(f"Mean normalized distance - with LCRs: {results['mean_normalized_distance_with_lcrs']:.3f}")
        print(f"Mean normalized distance - without LCRs: {results['mean_normalized_distance_without_lcrs']:.3f}")
        print(f"Mean genomic distance - with LCRs: {results['mean_genomic_distance_with_lcrs']:.0f} bp")
        print(f"Mean genomic distance - without LCRs: {results['mean_genomic_distance_without_lcrs']:.0f} bp")
        print(f"Mean gene length - with LCRs: {results['mean_gene_length_with_lcrs']:.0f} bp")
        print(f"Mean gene length - without LCRs: {results['mean_gene_length_without_lcrs']:.0f} bp")
        
        if distance_stats:
            print(f"Normalized distance comparison p-value: {distance_stats.pvalue:.3e}")
        if genomic_distance_stats:
            print(f"Genomic distance comparison p-value: {genomic_distance_stats.pvalue:.3e}")
        if length_stats:
            print(f"Gene length comparison p-value: {length_stats.pvalue:.3e}")
        if lcr_distance_correlation:
            print(f"LCR count vs distance correlation: r={lcr_distance_correlation[0]:.3f}, p={lcr_distance_correlation[1]:.3e}")
        
        return results#!/usr/bin/env python3
"""
LCR Paralog Analysis Script for SEG-identified Low Complexity Regions
Analyzes the corrected LCR data to answer three key questions:
1. LCR relative positions and overlap in paralog pairs
2. Evolutionary rates of genes with vs without LCRs
3. Evolutionary distances of LCRs and genes in paralogous pairs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

class LCRParalogAnalyzer:
    def __init__(self, data_file, paralog_coords_file=None, paralog_sequences_file=None):
        """
        Initialize analyzer with the corrected LCR data file and optional additional files
        
        Args:
            data_file: Path to corrected_lcr_data_with_chromosomes.v1.xlsb.csv
            paralog_coords_file: Path to paralog_coordinates.tsv (optional)
            paralog_sequences_file: Path to valid_paralog_nucleotide_sequences.fasta (optional)
        """
        self.data_file = data_file
        self.paralog_coords_file = paralog_coords_file
        self.paralog_sequences_file = paralog_sequences_file
        self.df = pd.read_csv(data_file)
        self.results = {}
        
        # Load additional data if provided
        self.paralog_coords_df = None
        self.all_genes_info = None
        self.sequence_genes = set()  # Genes from FASTA file
        
        if paralog_coords_file:
            try:
                self.paralog_coords_df = pd.read_csv(paralog_coords_file, sep='\t')
                print(f"Loaded paralog coordinates: {len(self.paralog_coords_df)} entries")
            except Exception as e:
                print(f"Could not load paralog coordinates file: {e}")
        
        if paralog_sequences_file:
            try:
                self.sequence_genes = self._load_sequence_gene_ids(paralog_sequences_file)
                print(f"Loaded sequence gene IDs: {len(self.sequence_genes)} genes")
            except Exception as e:
                print(f"Could not load sequences file: {e}")
        
        # Clean and prepare data
        self._prepare_data()
        self._load_all_genes_info()
    
    def _load_sequence_gene_ids(self, fasta_file):
        """Extract gene IDs from FASTA file headers"""
        gene_ids = set()
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Extract gene ID from header (remove '>' and take first part)
                        header = line.strip()[1:]  # Remove '>'
                        # Handle different header formats
                        gene_id = header.split()[0]  # Take first part before space
                        gene_id = gene_id.split('|')[0]  # Remove pipe separators if present
                        gene_ids.add(gene_id)
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
        return gene_ids
        
    def _prepare_data(self):
        """Clean and prepare the data for analysis"""
        print("Preparing data...")
        
        # Convert boolean columns if they're strings
        bool_cols = ['Same_Chromosome', 'Distance_Reliable']
        for col in bool_cols:
            if col in self.df.columns:
                if self.df[col].dtype == 'object':
                    self.df[col] = self.df[col].map({'True': True, 'FALSE': False, 'False': False, True: True, False: False})
        
        # Parse LCR positions to extract start and end positions
        self.df['Gene1_LCR_Start'] = np.nan
        self.df['Gene1_LCR_End'] = np.nan
        
        for idx, row in self.df.iterrows():
            if pd.notna(row['LCR_Positions_Corrected']) and row['LCR_Positions_Corrected'] != 'nan':
                positions = str(row['LCR_Positions_Corrected'])
                # Handle multiple LCR positions (take first one for simplicity)
                if ';' in positions:
                    positions = positions.split(';')[0]
                if '-' in positions:
                    try:
                        start, end = positions.split('-')
                        self.df.at[idx, 'Gene1_LCR_Start'] = int(start)
                        self.df.at[idx, 'Gene1_LCR_End'] = int(end)
                    except:
                        pass
        
        # Calculate relative LCR positions within genes
        self.df['Gene1_LCR_Relative_Start'] = (self.df['Gene1_LCR_Start'] - 1) / self.df['Gene1_Length_Genomic']  # -1 for 0-based
        self.df['Gene1_LCR_Relative_End'] = (self.df['Gene1_LCR_End'] - 1) / self.df['Gene1_Length_Genomic']
        self.df['Gene1_LCR_Relative_Mid'] = (self.df['Gene1_LCR_Relative_Start'] + self.df['Gene1_LCR_Relative_End']) / 2
        
        # Create binary indicators
        self.df['Has_LCR'] = self.df['Number_of_LCRs'] > 0
        self.df['Has_Paralog_Bool'] = self.df['Has_Paralog'] == 'Yes'
        self.df['Both_Have_LCR_Bool'] = self.df['Both_Have_LCR'] == 'Yes'
        
        # Create a lookup for paralog LCR positions
        self._create_paralog_lcr_positions()
        
        print(f"Data prepared: {len(self.df)} rows")
        print(f"Genes with LCRs: {self.df['Has_LCR'].sum()}")
        print(f"Genes with paralogs: {self.df['Has_Paralog_Bool'].sum()}")
        print(f"Paralog pairs where both have LCRs: {self.df['Both_Have_LCR_Bool'].sum()}")
    
    def _create_paralog_lcr_positions(self):
        """Create a lookup for paralog LCR positions"""
        # Create a dictionary to store LCR positions for each gene
        gene_lcr_positions = {}
        
        for idx, row in self.df.iterrows():
            if pd.notna(row['Gene1_LCR_Relative_Mid']):
                gene_lcr_positions[row['Gene_ID']] = {
                    'relative_start': row['Gene1_LCR_Relative_Start'],
                    'relative_end': row['Gene1_LCR_Relative_End'],
                    'relative_mid': row['Gene1_LCR_Relative_Mid'],
                    'length': row['LCR_Length_Corrected']
                }
        
        # Add paralog LCR positions
        self.df['Gene2_LCR_Relative_Start'] = np.nan
        self.df['Gene2_LCR_Relative_End'] = np.nan
        self.df['Gene2_LCR_Relative_Mid'] = np.nan
        self.df['Gene2_LCR_Length'] = np.nan
        
        for idx, row in self.df.iterrows():
            if pd.notna(row['Paralog_ID']) and row['Paralog_ID'] in gene_lcr_positions:
                paralog_info = gene_lcr_positions[row['Paralog_ID']]
                self.df.at[idx, 'Gene2_LCR_Relative_Start'] = paralog_info['relative_start']
                self.df.at[idx, 'Gene2_LCR_Relative_End'] = paralog_info['relative_end']
                self.df.at[idx, 'Gene2_LCR_Relative_Mid'] = paralog_info['relative_mid']
                self.df.at[idx, 'Gene2_LCR_Length'] = paralog_info['length']

    def _load_all_genes_info(self):
        """Load information about all genes from additional files"""
        self.all_genes_info = {}
        
        # First, get all genes from the LCR dataset (these all have LCRs)
        for idx, row in self.df.iterrows():
            # Gene1 info
            self.all_genes_info[row['Gene_ID']] = {
                'has_lcr': True,
                'lcr_count': row['Number_of_LCRs'],
                'gene_length': row['Gene1_Length_Genomic'],
                'chromosome': row['Gene1_Chromosome'],
                'start': row['Gene1_Start'],
                'end': row['Gene1_End'],
                'strand': row['Gene1_Strand']
            }
            
            # Gene2 (paralog) info if available
            if pd.notna(row['Paralog_ID']):
                if row['Paralog_ID'] not in self.all_genes_info:
                    self.all_genes_info[row['Paralog_ID']] = {
                        'has_lcr': row['Paralog_LCR_Count'] > 0,
                        'lcr_count': row['Paralog_LCR_Count'],
                        'gene_length': row['Gene2_Length_Genomic'],
                        'chromosome': row['Gene2_Chromosome'],
                        'start': row['Gene2_Start'],
                        'end': row['Gene2_End'],
                        'strand': row['Gene2_Strand']
                    }
        
        # Load genes from FASTA file that don't have LCRs
        if self.sequence_genes:
            print("Identifying genes without LCRs from sequence file...")
            genes_with_lcrs_set = set(gene_id for gene_id, info in self.all_genes_info.items() if info['has_lcr'])
            
            for gene_id in self.sequence_genes:
                if gene_id not in self.all_genes_info:
                    # This gene is in the FASTA but not in LCR dataset = no LCRs
                    # Try to get coordinates from paralog_coords if available
                    gene_length = 1000  # Default length
                    chromosome = 'Unknown'
                    start = 0
                    end = 1000
                    strand = '+'
                    
                    # Try to find coordinates in paralog_coords_df
                    if self.paralog_coords_df is not None:
                        # Find possible column names
                        possible_id_cols = ['gene_id', 'Gene_ID', 'ID', 'transcript_id', 'gene']
                        possible_chr_cols = ['chromosome', 'chr', 'chrom', 'Chromosome']
                        possible_start_cols = ['start', 'Start', 'gene_start', 'txStart']
                        possible_end_cols = ['end', 'End', 'gene_end', 'txEnd']
                        possible_strand_cols = ['strand', 'Strand']
                        
                        id_col = None
                        chr_col = None
                        start_col = None
                        end_col = None
                        strand_col = None
                        
                        for col in self.paralog_coords_df.columns:
                            if col in possible_id_cols and id_col is None:
                                id_col = col
                            elif col in possible_chr_cols and chr_col is None:
                                chr_col = col
                            elif col in possible_start_cols and start_col is None:
                                start_col = col
                            elif col in possible_end_cols and end_col is None:
                                end_col = col
                            elif col in possible_strand_cols and strand_col is None:
                                strand_col = col
                        
                        # Look for this gene in paralog coordinates
                        if id_col:
                            gene_row = self.paralog_coords_df[self.paralog_coords_df[id_col] == gene_id]
                            if not gene_row.empty:
                                row = gene_row.iloc[0]
                                if start_col and end_col and pd.notna(row[start_col]) and pd.notna(row[end_col]):
                                    start = int(row[start_col])
                                    end = int(row[end_col])
                                    gene_length = abs(end - start)
                                if chr_col and pd.notna(row[chr_col]):
                                    chromosome = row[chr_col]
                                if strand_col and pd.notna(row[strand_col]):
                                    strand = row[strand_col]
                    
                    self.all_genes_info[gene_id] = {
                        'has_lcr': False,
                        'lcr_count': 0,
                        'gene_length': gene_length,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }
        
        # Load additional genes from paralog coordinates file if available
        if self.paralog_coords_df is not None:
            print("Loading additional gene information from paralog coordinates...")
            
            # Try to identify column names (common variations)
            possible_id_cols = ['gene_id', 'Gene_ID', 'ID', 'transcript_id', 'gene']
            possible_chr_cols = ['chromosome', 'chr', 'chrom', 'Chromosome']
            possible_start_cols = ['start', 'Start', 'gene_start', 'txStart']
            possible_end_cols = ['end', 'End', 'gene_end', 'txEnd']
            possible_strand_cols = ['strand', 'Strand']
            
            # Find actual column names
            id_col = None
            chr_col = None
            start_col = None
            end_col = None
            strand_col = None
            
            for col in self.paralog_coords_df.columns:
                if col in possible_id_cols:
                    id_col = col
                elif col in possible_chr_cols:
                    chr_col = col
                elif col in possible_start_cols:
                    start_col = col
                elif col in possible_end_cols:
                    end_col = col
                elif col in possible_strand_cols:
                    strand_col = col
            
            print(f"Paralog coordinates columns: {list(self.paralog_coords_df.columns)}")
            print(f"Using: ID={id_col}, Chr={chr_col}, Start={start_col}, End={end_col}, Strand={strand_col}")
            
            # Add genes from paralog coordinates that don't have LCRs
            if id_col:
                for idx, row in self.paralog_coords_df.iterrows():
                    gene_id = row[id_col]
                    if gene_id not in self.all_genes_info:
                        # This gene doesn't have LCRs (not in our LCR dataset)
                        gene_length = 0
                        if start_col and end_col and pd.notna(row[start_col]) and pd.notna(row[end_col]):
                            gene_length = abs(int(row[end_col]) - int(row[start_col]))
                        
                        self.all_genes_info[gene_id] = {
                            'has_lcr': False,
                            'lcr_count': 0,
                            'gene_length': gene_length,
                            'chromosome': row[chr_col] if chr_col and pd.notna(row[chr_col]) else 'Unknown',
                            'start': int(row[start_col]) if start_col and pd.notna(row[start_col]) else 0,
                            'end': int(row[end_col]) if end_col and pd.notna(row[end_col]) else 0,
                            'strand': row[strand_col] if strand_col and pd.notna(row[strand_col]) else '+'
                        }
        
        genes_with_lcrs = sum(1 for g in self.all_genes_info.values() if g['has_lcr'])
        genes_without_lcrs = sum(1 for g in self.all_genes_info.values() if not g['has_lcr'])
        
        print(f"Total genes loaded: {len(self.all_genes_info)}")
        print(f"Genes with LCRs: {genes_with_lcrs}")
        print(f"Genes without LCRs: {genes_without_lcrs}")
        
        return self.all_genes_info
    
    def question1_lcr_positions_overlap(self):
        """
        Question 1: Calculate relative positions of LCRs in paralog pairs
        to compare LCR overlap and distinguish related vs independent formation
        """
        print("\n=== QUESTION 1: LCR Positions and Overlap Analysis ===")
        
        # Filter for paralog pairs where both have LCRs and we have position data
        paralog_pairs = self.df[
            (self.df['Has_Paralog_Bool'] == True) & 
            (self.df['Both_Have_LCR_Bool'] == True) &
            (pd.notna(self.df['Gene1_LCR_Relative_Mid'])) &
            (pd.notna(self.df['Gene2_LCR_Relative_Mid']))
        ].copy()
        
        print(f"Analyzing {len(paralog_pairs)} paralog pairs with LCR position data")
        
        analysis_results = []
        
        for idx, row in paralog_pairs.iterrows():
            # LCR relative positions within each gene
            gene1_lcr_rel_start = row['Gene1_LCR_Relative_Start']
            gene1_lcr_rel_end = row['Gene1_LCR_Relative_End']
            gene1_lcr_rel_mid = row['Gene1_LCR_Relative_Mid']
            
            gene2_lcr_rel_start = row['Gene2_LCR_Relative_Start']
            gene2_lcr_rel_end = row['Gene2_LCR_Relative_End']
            gene2_lcr_rel_mid = row['Gene2_LCR_Relative_Mid']
            
            # Calculate position similarities
            start_position_similarity = 1 - abs(gene1_lcr_rel_start - gene2_lcr_rel_start)
            end_position_similarity = 1 - abs(gene1_lcr_rel_end - gene2_lcr_rel_end)
            mid_position_similarity = 1 - abs(gene1_lcr_rel_mid - gene2_lcr_rel_mid)
            
            # Overall position similarity (average of start, end, mid)
            overall_position_similarity = (start_position_similarity + end_position_similarity + mid_position_similarity) / 3
            
            # LCR length similarity
            gene1_lcr_length = row['LCR_Length_Corrected']
            gene2_lcr_length = row['Gene2_LCR_Length'] if pd.notna(row['Gene2_LCR_Length']) else gene1_lcr_length
            lcr_length_similarity = 1 - abs(gene1_lcr_length - gene2_lcr_length) / max(gene1_lcr_length, gene2_lcr_length, 1)
            
            # Calculate LCR count similarity
            lcr_count_similarity = 1 - abs(row['Number_of_LCRs'] - row['Paralog_LCR_Count']) / max(row['Number_of_LCRs'], row['Paralog_LCR_Count'], 1)
            
            # Physical overlap potential (only for same chromosome)
            physical_overlap = 0
            gene1_lcr_genomic_start = 0
            gene1_lcr_genomic_end = 0
            gene2_lcr_genomic_start = 0
            gene2_lcr_genomic_end = 0
            
            if row['Same_Chromosome']:
                # Convert relative positions back to genomic coordinates
                gene1_lcr_genomic_start = row['Gene1_Start'] + (gene1_lcr_rel_start * row['Gene1_Length_Genomic'])
                gene1_lcr_genomic_end = row['Gene1_Start'] + (gene1_lcr_rel_end * row['Gene1_Length_Genomic'])
                gene2_lcr_genomic_start = row['Gene2_Start'] + (gene2_lcr_rel_start * row['Gene2_Length_Genomic'])
                gene2_lcr_genomic_end = row['Gene2_Start'] + (gene2_lcr_rel_end * row['Gene2_Length_Genomic'])
                
                # Calculate actual overlap
                overlap_start = max(gene1_lcr_genomic_start, gene2_lcr_genomic_start)
                overlap_end = min(gene1_lcr_genomic_end, gene2_lcr_genomic_end)
                physical_overlap = max(0, overlap_end - overlap_start)
            else:
                # Still calculate genomic coordinates even if different chromosomes
                gene1_lcr_genomic_start = row['Gene1_Start'] + (gene1_lcr_rel_start * row['Gene1_Length_Genomic'])
                gene1_lcr_genomic_end = row['Gene1_Start'] + (gene1_lcr_rel_end * row['Gene1_Length_Genomic'])
                gene2_lcr_genomic_start = row['Gene2_Start'] + (gene2_lcr_rel_start * row['Gene2_Length_Genomic'])
                gene2_lcr_genomic_end = row['Gene2_Start'] + (gene2_lcr_rel_end * row['Gene2_Length_Genomic'])
            
            # Determine if LCRs likely formed together (related) or independently
            # Criteria: high position similarity, same chromosome, similar lengths
            is_related_formation = (
                (overall_position_similarity > 0.7) and 
                (lcr_length_similarity > 0.5) and
                (lcr_count_similarity > 0.7) and
                ((row['Same_Chromosome'] and row['Genomic_Distance_bp'] < 1000000) or overall_position_similarity > 0.9)
            )
            
            analysis_results.append({
                'Gene1_ID': row['Gene_ID'],
                'Gene2_ID': row['Paralog_ID'],
                'Gene1_Chromosome': row['Gene1_Chromosome'],
                'Gene2_Chromosome': row['Gene2_Chromosome'],
                'Same_Chromosome': row['Same_Chromosome'],
                'Genomic_Distance_bp': row['Genomic_Distance_bp'],
                
                # Gene coordinates
                'Gene1_Start': row['Gene1_Start'],
                'Gene1_End': row['Gene1_End'],
                'Gene2_Start': row['Gene2_Start'],
                'Gene2_End': row['Gene2_End'],
                
                # LCR relative positions
                'Gene1_LCR_Relative_Start': gene1_lcr_rel_start,
                'Gene1_LCR_Relative_End': gene1_lcr_rel_end,
                'Gene1_LCR_Relative_Mid': gene1_lcr_rel_mid,
                'Gene2_LCR_Relative_Start': gene2_lcr_rel_start,
                'Gene2_LCR_Relative_End': gene2_lcr_rel_end,
                'Gene2_LCR_Relative_Mid': gene2_lcr_rel_mid,
                
                # LCR genomic coordinates
                'Gene1_LCR_Genomic_Start': int(gene1_lcr_genomic_start),
                'Gene1_LCR_Genomic_End': int(gene1_lcr_genomic_end),
                'Gene2_LCR_Genomic_Start': int(gene2_lcr_genomic_start),
                'Gene2_LCR_Genomic_End': int(gene2_lcr_genomic_end),
                
                # LCR lengths
                'Gene1_LCR_Length': gene1_lcr_length,
                'Gene2_LCR_Length': gene2_lcr_length,
                
                # Similarity metrics
                'Start_Position_Similarity': start_position_similarity,
                'End_Position_Similarity': end_position_similarity,
                'Mid_Position_Similarity': mid_position_similarity,
                'Overall_Position_Similarity': overall_position_similarity,
                'LCR_Length_Similarity': lcr_length_similarity,
                'LCR_Count_Similarity': lcr_count_similarity,
                'Physical_Overlap_bp': physical_overlap,
                'Is_Related_Formation': is_related_formation,
                'Distance_Category': row['Distance_Category']
            })
        
        results_df = pd.DataFrame(analysis_results)
        
        # Summary statistics
        if len(results_df) > 0:
            related_count = results_df['Is_Related_Formation'].sum()
            independent_count = len(results_df) - related_count
            
            summary_stats = {
                'total_paralog_pairs': len(results_df),
                'related_formations': related_count,
                'independent_formations': independent_count,
                'percent_related': (related_count / len(results_df)) * 100,
                'mean_overall_position_similarity': results_df['Overall_Position_Similarity'].mean(),
                'mean_lcr_length_similarity': results_df['LCR_Length_Similarity'].mean(),
                'mean_physical_overlap': results_df['Physical_Overlap_bp'].mean(),
                'same_chromosome_pairs': results_df['Same_Chromosome'].sum(),
                'mean_genomic_distance': results_df['Genomic_Distance_bp'].mean()
            }
        else:
            summary_stats = {
                'total_paralog_pairs': 0,
                'related_formations': 0,
                'independent_formations': 0,
                'percent_related': 0,
                'mean_overall_position_similarity': 0,
                'mean_lcr_length_similarity': 0,
                'mean_physical_overlap': 0,
                'same_chromosome_pairs': 0,
                'mean_genomic_distance': 0
            }
        
        self.results['question1'] = {
            'detailed_results': results_df,
            'summary_stats': summary_stats
        }
        
        print(f"Related LCR formations: {summary_stats['related_formations']} ({summary_stats['percent_related']:.1f}%)")
        print(f"Independent LCR formations: {summary_stats['independent_formations']}")
        print(f"Mean overall position similarity: {summary_stats['mean_overall_position_similarity']:.3f}")
        print(f"Mean physical overlap: {summary_stats['mean_physical_overlap']:.1f} bp")
        
        return results_df, summary_stats
    
    def question2_genetic_distances(self):
        """
        Question 2: Calculate genetic distances of genes with vs without LCRs
        to compare their evolutionary rates
        """
        print("\n=== QUESTION 2: Genetic Distance Analysis ===")
        
        # Use comprehensive gene information
        if self.all_genes_info is None:
            print("No comprehensive gene information available. Using only LCR dataset.")
            genes_with_lcrs = [g for g, info in self.all_genes_info.items() if info['has_lcr']] if self.all_genes_info else []
            genes_without_lcrs = []
        else:
            genes_with_lcrs = [g for g, info in self.all_genes_info.items() if info['has_lcr']]
            genes_without_lcrs = [g for g, info in self.all_genes_info.items() if not info['has_lcr']]
        
        print(f"Total unique genes: {len(self.all_genes_info) if self.all_genes_info else 0}")
        print(f"Genes with LCRs: {len(genes_with_lcrs)}")
        print(f"Genes without LCRs: {len(genes_without_lcrs)}")
        
        # Calculate evolutionary distances using paralog pairs from the LCR dataset
        def calculate_distances_for_group(gene_list, group_name):
            distances = []
            for idx, row in self.df.iterrows():
                if (row['Gene_ID'] in gene_list and 
                    row['Has_Paralog_Bool'] and 
                    row['Distance_Reliable'] and
                    pd.notna(row['Genomic_Distance_bp'])):
                    
                    # Use genomic distance normalized by gene length as evolutionary proxy
                    normalized_distance = row['Genomic_Distance_bp'] / row['Gene1_Length_Genomic'] if row['Gene1_Length_Genomic'] > 0 else 0
                    distances.append({
                        'gene_id': row['Gene_ID'],
                        'paralog_id': row['Paralog_ID'],
                        'genomic_distance': row['Genomic_Distance_bp'],
                        'normalized_distance': normalized_distance,
                        'same_chromosome': row['Same_Chromosome'],
                        'gene_length': row['Gene1_Length_Genomic'],
                        'lcr_count': row['Number_of_LCRs']
                    })
            print(f"Found {len(distances)} paralog pairs for genes {group_name}")
            return distances
        
        # Calculate distances for both groups
        lcr_distances = calculate_distances_for_group(genes_with_lcrs, "with LCRs")
        no_lcr_distances = calculate_distances_for_group(genes_without_lcrs, "without LCRs")
        
        # Extract numerical values for statistical analysis
        lcr_distance_values = [d['normalized_distance'] for d in lcr_distances if d['normalized_distance'] > 0]
        no_lcr_distance_values = [d['normalized_distance'] for d in no_lcr_distances if d['normalized_distance'] > 0]
        
        lcr_genomic_values = [d['genomic_distance'] for d in lcr_distances]
        no_lcr_genomic_values = [d['genomic_distance'] for d in no_lcr_distances]
        
        # Calculate gene length statistics using comprehensive gene info
        lcr_gene_lengths = []
        no_lcr_gene_lengths = []
        
        if self.all_genes_info:
            lcr_gene_lengths = [info['gene_length'] for gene_id, info in self.all_genes_info.items() 
                              if info['has_lcr'] and info['gene_length'] > 0]
            no_lcr_gene_lengths = [info['gene_length'] for gene_id, info in self.all_genes_info.items() 
                                 if not info['has_lcr'] and info['gene_length'] > 0]
        
        # Statistical comparisons
        distance_stats = None
        genomic_distance_stats = None
        length_stats = None
        
        if len(lcr_distance_values) > 0 and len(no_lcr_distance_values) > 0:
            distance_stats = mannwhitneyu(lcr_distance_values, no_lcr_distance_values, alternative='two-sided')
            
        if len(lcr_genomic_values) > 0 and len(no_lcr_genomic_values) > 0:
            genomic_distance_stats = mannwhitneyu(lcr_genomic_values, no_lcr_genomic_values, alternative='two-sided')
        
        if len(lcr_gene_lengths) > 0 and len(no_lcr_gene_lengths) > 0:
            length_stats = mannwhitneyu(lcr_gene_lengths, no_lcr_gene_lengths, alternative='two-sided')
        
        # Correlation analysis: LCR count vs evolutionary distance
        lcr_distance_correlation = None
        if len(lcr_distances) > 1:
            lcr_counts = [d['lcr_count'] for d in lcr_distances]
            distances = [d['normalized_distance'] for d in lcr_distances]
            # Filter out zero distances for correlation
            valid_pairs = [(lc, dist) for lc, dist in zip(lcr_counts, distances) if dist > 0]
            if len(valid_pairs) > 1:
                lcr_counts_valid = [x[0] for x in valid_pairs]
                distances_valid = [x[1] for x in valid_pairs]
                lcr_distance_correlation = pearsonr(lcr_counts_valid, distances_valid)
        
        results = {
            'total_genes_analyzed': len(self.all_genes_info) if self.all_genes_info else len(set(self.df['Gene_ID'])),
            'genes_with_lcrs_count': len(genes_with_lcrs),
            'genes_without_lcrs_count': len(genes_without_lcrs),
            'lcr_paralog_pairs_analyzed': len(lcr_distances),
            'no_lcr_paralog_pairs_analyzed': len(no_lcr_distances),
            'mean_normalized_distance_with_lcrs': np.mean(lcr_distance_values) if lcr_distance_values else 0,
            'mean_normalized_distance_without_lcrs': np.mean(no_lcr_distance_values) if no_lcr_distance_values else 0,
            'median_normalized_distance_with_lcrs': np.median(lcr_distance_values) if lcr_distance_values else 0,
            'median_normalized_distance_without_lcrs': np.median(no_lcr_distance_values) if no_lcr_distance_values else 0,
            'mean_genomic_distance_with_lcrs': np.mean(lcr_genomic_values) if lcr_genomic_values else 0,
            'mean_genomic_distance_without_lcrs': np.mean(no_lcr_genomic_values) if no_lcr_genomic_values else 0,
            'mean_gene_length_with_lcrs': np.mean(lcr_gene_lengths) if lcr_gene_lengths else 0,
            'mean_gene_length_without_lcrs': np.mean(no_lcr_gene_lengths) if no_lcr_gene_lengths else 0,
            'normalized_distance_statistical_test': distance_stats,
            'genomic_distance_statistical_test': genomic_distance_stats,
            'gene_length_statistical_test': length_stats,
            'lcr_count_distance_correlation': lcr_distance_correlation,
            'detailed_lcr_distances': lcr_distances,
            'detailed_no_lcr_distances': no_lcr_distances
        }
        
        self.results['question2'] = results
        
        print(f"Paralog pairs analyzed - with LCRs: {len(lcr_distances)}")
        print(f"Paralog pairs analyzed - without LCRs: {len(no_lcr_distances)}")
        print(f"Mean normalized distance - with LCRs: {results['mean_normalized_distance_with_lcrs']:.3f}")
        print(f"Mean normalized distance - without LCRs: {results['mean_normalized_distance_without_lcrs']:.3f}")
        print(f"Mean genomic distance - with LCRs: {results['mean_genomic_distance_with_lcrs']:.0f} bp")
        print(f"Mean genomic distance - without LCRs: {results['mean_genomic_distance_without_lcrs']:.0f} bp")
        print(f"Mean gene length - with LCRs: {results['mean_gene_length_with_lcrs']:.0f} bp")
        print(f"Mean gene length - without LCRs: {results['mean_gene_length_without_lcrs']:.0f} bp")
        
        if distance_stats:
            print(f"Normalized distance comparison p-value: {distance_stats.pvalue:.3e}")
        if genomic_distance_stats:
            print(f"Genomic distance comparison p-value: {genomic_distance_stats.pvalue:.3e}")
        if length_stats:
            print(f"Gene length comparison p-value: {length_stats.pvalue:.3e}")
        if lcr_distance_correlation:
            print(f"LCR count vs distance correlation: r={lcr_distance_correlation[0]:.3f}, p={lcr_distance_correlation[1]:.3e}")
        
        return results
    
    def question3_evolutionary_distances(self):
        """
        Question 3: Calculate evolutionary distances of LCRs and genes in paralogous pairs
        to compare with recent duplications
        """
        print("\n=== QUESTION 3: Evolutionary Distance Analysis ===")
        
        # Focus on paralog pairs with reliable distance measurements
        paralog_pairs = self.df[
            (self.df['Has_Paralog_Bool'] == True) & 
            (self.df['Distance_Reliable'] == True)
        ].copy()
        
        print(f"Analyzing {len(paralog_pairs)} paralog pairs with reliable distance data")
        
        # Classify duplications by distance and chromosome location
        paralog_pairs['Duplication_Type'] = 'Ancient'
        
        # Recent duplications: same chromosome, close distance
        recent_mask = (
            (paralog_pairs['Same_Chromosome'] == True) & 
            (paralog_pairs['Genomic_Distance_bp'] < 100000)  # Within 100kb
        )
        paralog_pairs.loc[recent_mask, 'Duplication_Type'] = 'Recent'
        
        # Intermediate duplications: same chromosome, moderate distance
        intermediate_mask = (
            (paralog_pairs['Same_Chromosome'] == True) & 
            (paralog_pairs['Genomic_Distance_bp'] >= 100000) & 
            (paralog_pairs['Genomic_Distance_bp'] < 1000000)  # 100kb - 1Mb
        )
        paralog_pairs.loc[intermediate_mask, 'Duplication_Type'] = 'Intermediate'
        
        # Calculate evolutionary distance proxies
        evolutionary_analysis = []
        
        for idx, row in paralog_pairs.iterrows():
            # Gene evolutionary distance proxy (normalized genomic distance)
            gene_distance = row['Genomic_Distance_bp'] / max(row['Gene1_Length_Genomic'], row['Gene2_Length_Genomic'])
            
            # LCR evolutionary distance proxy
            # Based on difference in LCR counts and patterns
            lcr_count_diff = abs(row['Number_of_LCRs'] - row['Paralog_LCR_Count'])
            max_lcr_count = max(row['Number_of_LCRs'], row['Paralog_LCR_Count'])
            lcr_divergence = lcr_count_diff / max_lcr_count if max_lcr_count > 0 else 0
            
            # LCR vs gene distance ratio
            lcr_gene_ratio = lcr_divergence / gene_distance if gene_distance > 0 else 0
            
            evolutionary_analysis.append({
                'Gene1_ID': row['Gene_ID'],
                'Gene2_ID': row['Paralog_ID'],
                'Duplication_Type': row['Duplication_Type'],
                'Genomic_Distance_bp': row['Genomic_Distance_bp'],
                'Gene_Distance_Normalized': gene_distance,
                'LCR_Divergence': lcr_divergence,
                'LCR_Gene_Ratio': lcr_gene_ratio,
                'Gene1_LCR_Count': row['Number_of_LCRs'],
                'Gene2_LCR_Count': row['Paralog_LCR_Count'],
                'Same_Chromosome': row['Same_Chromosome'],
                'Both_Have_LCR': row['Both_Have_LCR_Bool']
            })
        
        evo_df = pd.DataFrame(evolutionary_analysis)
        
        # Compare evolutionary distances by duplication type
        recent_pairs = evo_df[evo_df['Duplication_Type'] == 'Recent']
        intermediate_pairs = evo_df[evo_df['Duplication_Type'] == 'Intermediate']
        ancient_pairs = evo_df[evo_df['Duplication_Type'] == 'Ancient']
        
        # Statistical comparisons
        comparisons = {}
        
        # Recent vs Ancient
        if len(recent_pairs) > 0 and len(ancient_pairs) > 0:
            comparisons['recent_vs_ancient_gene'] = mannwhitneyu(
                recent_pairs['Gene_Distance_Normalized'], 
                ancient_pairs['Gene_Distance_Normalized']
            )
            comparisons['recent_vs_ancient_lcr'] = mannwhitneyu(
                recent_pairs['LCR_Divergence'], 
                ancient_pairs['LCR_Divergence']
            )
        
        # LCR vs Gene evolution rates
        lcr_gene_correlation = pearsonr(evo_df['Gene_Distance_Normalized'], evo_df['LCR_Divergence'])
        
        summary_stats = {
            'total_pairs': len(evo_df),
            'recent_duplications': len(recent_pairs),
            'intermediate_duplications': len(intermediate_pairs),
            'ancient_duplications': len(ancient_pairs),
            'mean_gene_distance_recent': recent_pairs['Gene_Distance_Normalized'].mean() if len(recent_pairs) > 0 else 0,
            'mean_gene_distance_ancient': ancient_pairs['Gene_Distance_Normalized'].mean() if len(ancient_pairs) > 0 else 0,
            'mean_lcr_divergence_recent': recent_pairs['LCR_Divergence'].mean() if len(recent_pairs) > 0 else 0,
            'mean_lcr_divergence_ancient': ancient_pairs['LCR_Divergence'].mean() if len(ancient_pairs) > 0 else 0,
            'lcr_gene_correlation': lcr_gene_correlation,
            'statistical_comparisons': comparisons
        }
        
        self.results['question3'] = {
            'detailed_results': evo_df,
            'summary_stats': summary_stats
        }
        
        print(f"Recent duplications: {len(recent_pairs)}")
        print(f"Intermediate duplications: {len(intermediate_pairs)}")
        print(f"Ancient duplications: {len(ancient_pairs)}")
        print(f"Mean gene distance (recent): {summary_stats['mean_gene_distance_recent']:.3f}")
        print(f"Mean gene distance (ancient): {summary_stats['mean_gene_distance_ancient']:.3f}")
        print(f"Mean LCR divergence (recent): {summary_stats['mean_lcr_divergence_recent']:.3f}")
        print(f"Mean LCR divergence (ancient): {summary_stats['mean_lcr_divergence_ancient']:.3f}")
        print(f"LCR-Gene correlation: r={lcr_gene_correlation[0]:.3f}, p={lcr_gene_correlation[1]:.3e}")
        
        return evo_df, summary_stats
    
    def generate_comprehensive_summary(self):
        """Generate a comprehensive summary answering all three questions"""
        print("\n=== COMPREHENSIVE SUMMARY ===")
        
        summary = {}
        
        # Question 1 Summary
        if 'question1' in self.results:
            q1_stats = self.results['question1']['summary_stats']
            summary['Q1_LCR_Formation'] = {
                'total_paralog_pairs_analyzed': q1_stats['total_paralog_pairs'],
                'related_formations_count': q1_stats['related_formations'],
                'independent_formations_count': q1_stats['independent_formations'],
                'percent_related_formation': q1_stats['percent_related'],
                'mean_overall_position_similarity': q1_stats['mean_overall_position_similarity'],
                'mean_lcr_length_similarity': q1_stats['mean_lcr_length_similarity'],
                'mean_physical_overlap': q1_stats['mean_physical_overlap'],
                'interpretation': 'Related formation' if q1_stats['percent_related'] > 50 else 'Independent formation'
            }
        
        # Question 2 Summary
        if 'question2' in self.results:
            q2_stats = self.results['question2']
            summary['Q2_Evolutionary_Rates'] = {
                'total_genes_analyzed': q2_stats['total_genes_analyzed'],
                'genes_with_lcrs': q2_stats['genes_with_lcrs_count'],
                'genes_without_lcrs': q2_stats['genes_without_lcrs_count'],
                'lcr_paralog_pairs_analyzed': q2_stats['lcr_paralog_pairs_analyzed'],
                'no_lcr_paralog_pairs_analyzed': q2_stats['no_lcr_paralog_pairs_analyzed'],
                'mean_normalized_distance_with_lcrs': q2_stats['mean_normalized_distance_with_lcrs'],
                'mean_normalized_distance_without_lcrs': q2_stats['mean_normalized_distance_without_lcrs'],
                'mean_genomic_distance_with_lcrs': q2_stats['mean_genomic_distance_with_lcrs'],
                'mean_genomic_distance_without_lcrs': q2_stats['mean_genomic_distance_without_lcrs'],
                'distance_difference': q2_stats['mean_normalized_distance_with_lcrs'] - q2_stats['mean_normalized_distance_without_lcrs'],
                'faster_evolution': 'With LCRs' if q2_stats['mean_normalized_distance_with_lcrs'] > q2_stats['mean_normalized_distance_without_lcrs'] else 'Without LCRs'
            }
        
        # Question 3 Summary
        if 'question3' in self.results:
            q3_stats = self.results['question3']['summary_stats']
            summary['Q3_Evolutionary_Distances'] = {
                'total_pairs_analyzed': q3_stats['total_pairs'],
                'recent_duplications': q3_stats['recent_duplications'],
                'ancient_duplications': q3_stats['ancient_duplications'],
                'recent_gene_distance': q3_stats['mean_gene_distance_recent'],
                'ancient_gene_distance': q3_stats['mean_gene_distance_ancient'],
                'recent_lcr_divergence': q3_stats['mean_lcr_divergence_recent'],
                'ancient_lcr_divergence': q3_stats['mean_lcr_divergence_ancient'],
                'lcr_gene_correlation': q3_stats['lcr_gene_correlation'][0] if q3_stats['lcr_gene_correlation'] else 0
            }
        
        return summary
    
    def save_results(self, output_prefix='lcr_analysis_results'):
        """Save all results to CSV files"""
        
        # Save detailed results for each question
        if 'question1' in self.results:
            self.results['question1']['detailed_results'].to_csv(
                f'{output_prefix}_Q1_lcr_positions.csv', index=False
            )
        
        if 'question3' in self.results:
            self.results['question3']['detailed_results'].to_csv(
                f'{output_prefix}_Q3_evolutionary_distances.csv', index=False
            )
        
        # Save question 2 detailed results
        if 'question2' in self.results:
            q2_results = self.results['question2']
            
            # Save detailed distance data
            if q2_results['detailed_lcr_distances']:
                lcr_dist_df = pd.DataFrame(q2_results['detailed_lcr_distances'])
                lcr_dist_df['gene_group'] = 'with_LCRs'
                lcr_dist_df.to_csv(f'{output_prefix}_Q2_detailed_lcr_distances.csv', index=False)
            
            if q2_results['detailed_no_lcr_distances']:
                no_lcr_dist_df = pd.DataFrame(q2_results['detailed_no_lcr_distances'])
                no_lcr_dist_df['gene_group'] = 'without_LCRs'
                no_lcr_dist_df.to_csv(f'{output_prefix}_Q2_detailed_no_lcr_distances.csv', index=False)
            
            # Save summary statistics
            q2_summary = {k: v for k, v in q2_results.items() 
                         if k not in ['detailed_lcr_distances', 'detailed_no_lcr_distances']}
            q2_df = pd.DataFrame([q2_summary])
            q2_df.to_csv(f'{output_prefix}_Q2_genetic_distances_summary.csv', index=False)
        
        # Save comprehensive summary
        summary = self.generate_comprehensive_summary()
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(f'{output_prefix}_comprehensive_summary.csv', index=False)
        
        print(f"\nResults saved with prefix: {output_prefix}")
        print("Files created:")
        print(f"- {output_prefix}_Q1_lcr_positions.csv")
        print(f"- {output_prefix}_Q2_detailed_lcr_distances.csv")
        print(f"- {output_prefix}_Q2_detailed_no_lcr_distances.csv") 
        print(f"- {output_prefix}_Q2_genetic_distances_summary.csv")
        print(f"- {output_prefix}_Q3_evolutionary_distances.csv")
        print(f"- {output_prefix}_comprehensive_summary.csv")
        
    def run_complete_analysis(self):
        """Run the complete analysis for all three questions"""
        print("Starting Complete LCR Paralog Analysis...")
        print(f"Dataset: {len(self.df)} total entries")
        
        # Run all three analyses
        q1_results = self.question1_lcr_positions_overlap()
        q2_results = self.question2_genetic_distances()
        q3_results = self.question3_evolutionary_distances()
        
        # Generate comprehensive summary
        summary = self.generate_comprehensive_summary()
        
        # Save results
        self.save_results()
        
        return summary

def main():
    """Main function to run the analysis"""
    import argparse
    
    parser = argparse.ArgumentParser(description='LCR Paralog Analysis for SEG-identified regions')
    parser.add_argument('--input', default='corrected_lcr_data_with_chromosomes.v1.xlsb.csv', 
                       help='Input CSV file with LCR data')
    parser.add_argument('--paralog-coords', default=None,
                       help='TSV file with paralog coordinates (paralog_coordinates.tsv)')
    parser.add_argument('--paralog-sequences', default=None,
                       help='FASTA file with paralog sequences (valid_paralog_nucleotide_sequences.fasta)')
    parser.add_argument('--output', default='lcr_analysis_results', 
                       help='Output prefix for result files')
    
    args = parser.parse_args()
    
    # Run analysis
    try:
        analyzer = LCRParalogAnalyzer(args.input, args.paralog_coords, args.paralog_sequences)
        summary = analyzer.run_complete_analysis()
        
        print("\n" + "="*50)
        print("ANALYSIS COMPLETE")
        print("="*50)
        print("Files generated:")
        print(f"- {args.output}_Q1_lcr_positions.csv")
        print(f"- {args.output}_Q2_detailed_lcr_distances.csv")
        print(f"- {args.output}_Q2_detailed_no_lcr_distances.csv") 
        print(f"- {args.output}_Q2_genetic_distances_summary.csv")
        print(f"- {args.output}_Q3_evolutionary_distances.csv")
        print(f"- {args.output}_comprehensive_summary.csv")
        
    except Exception as e:
        print(f"Error running analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
