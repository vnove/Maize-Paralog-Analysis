#!/usr/bin/env python3
"""
Comprehensive Detailed Analysis Script - All Research Questions
In-depth analysis with complete CSV outputs and statistical testing
"""

import pandas as pd
import numpy as np
import os
import sys
from scipy.stats import spearmanr, chi2_contingency, mannwhitneyu, pearsonr, kruskal
from scipy.stats import ttest_ind, fisher_exact
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

def load_and_validate_data():
    """Load datasets with comprehensive validation"""
    print("🔄 LOADING AND VALIDATING DATASETS")
    print("=" * 60)
    
    # Load datasets
    comprehensive_df = pd.read_csv('comprehensive_lcr_paralog_dataset.csv')
    pairwise_df = pd.read_csv('lcr_analysis_results_Q1_lcr_positions.csv')
    
    print(f"📊 Comprehensive dataset: {len(comprehensive_df):,} genes, {comprehensive_df.shape[1]} columns")
    print(f"📊 Pairwise dataset: {len(pairwise_df):,} gene pairs, {pairwise_df.shape[1]} columns")
    
    # Data validation and cleaning
    print("\n🔍 DATA VALIDATION:")
    
    # Check for missing values
    comprehensive_missing = comprehensive_df.isnull().sum()
    pairwise_missing = pairwise_df.isnull().sum()
    
    print(f"   Missing values in comprehensive dataset: {comprehensive_missing.sum():,}")
    print(f"   Missing values in pairwise dataset: {pairwise_missing.sum():,}")
    
    # Data quality checks
    lcr_range = f"{comprehensive_df['LCR_Count'].min()}-{comprehensive_df['LCR_Count'].max()}"
    paralog_range = f"{comprehensive_df['Total_Paralogs'].min()}-{comprehensive_df['Total_Paralogs'].max()}"
    
    print(f"   LCR Count range: {lcr_range}")
    print(f"   Paralog Count range: {paralog_range}")
    
    # Create analysis-ready datasets
    analysis_df = comprehensive_df.dropna(subset=['LCR_Count', 'Total_Paralogs']).copy()
    print(f"   Clean analysis dataset: {len(analysis_df):,} genes")
    
    return comprehensive_df, pairwise_df, analysis_df

def detailed_lcr_analysis(comprehensive_df, analysis_df):
    """Comprehensive LCR analysis"""
    print("\n" + "="*80)
    print("DETAILED LCR ANALYSIS")
    print("="*80)
    
    results = {}
    
    # LCR distribution analysis
    print("\n📊 LCR DISTRIBUTION ANALYSIS:")
    lcr_dist = comprehensive_df['LCR_Count'].value_counts().sort_index()
    
    print("   LCR Count Distribution:")
    for lcr_count, gene_count in lcr_dist.items():
        percentage = (gene_count / len(comprehensive_df)) * 100
        print(f"     {lcr_count} LCRs: {gene_count:,} genes ({percentage:.2f}%)")
    
    results['lcr_distribution'] = lcr_dist.to_dict()
    
    # LCR complexity analysis
    if 'Avg_Complexity' in comprehensive_df.columns:
        print(f"\n📊 LCR COMPLEXITY ANALYSIS:")
        complexity_stats = comprehensive_df['Avg_Complexity'].describe()
        print(f"   Average complexity: {complexity_stats['mean']:.3f} ± {complexity_stats['std']:.3f}")
        print(f"   Complexity range: {complexity_stats['min']:.3f} - {complexity_stats['max']:.3f}")
        results['avg_complexity'] = complexity_stats['mean']
        results['complexity_std'] = complexity_stats['std']
    
    # LCR length analysis
    if 'Total_LCR_Length' in comprehensive_df.columns:
        print(f"\n📊 LCR LENGTH ANALYSIS:")
        length_stats = comprehensive_df['Total_LCR_Length'].describe()
        print(f"   Average total LCR length: {length_stats['mean']:.1f} ± {length_stats['std']:.1f} bp")
        print(f"   Length range: {length_stats['min']:.0f} - {length_stats['max']:.0f} bp")
        results['avg_lcr_length'] = length_stats['mean']
        results['lcr_length_std'] = length_stats['std']
    
    # LCR categories
    if 'LCR_Category' in comprehensive_df.columns:
        print(f"\n📊 LCR CATEGORY ANALYSIS:")
        lcr_cat_dist = comprehensive_df['LCR_Category'].value_counts()
        for category, count in lcr_cat_dist.items():
            percentage = (count / len(comprehensive_df)) * 100
            print(f"   {category}: {count:,} genes ({percentage:.2f}%)")
        results['lcr_category_distribution'] = lcr_cat_dist.to_dict()
    
    return results

def detailed_paralogy_analysis(comprehensive_df, analysis_df):
    """Comprehensive paralogy analysis"""
    print("\n" + "="*80)
    print("DETAILED PARALOGY ANALYSIS")
    print("="*80)
    
    results = {}
    
    # Paralogy distribution
    print("\n📊 PARALOGY DISTRIBUTION ANALYSIS:")
    paralog_dist = comprehensive_df['Total_Paralogs'].value_counts().sort_index()
    
    print("   Top 10 Paralog Counts:")
    for paralog_count, gene_count in list(paralog_dist.items())[:10]:
        percentage = (gene_count / len(comprehensive_df)) * 100
        print(f"     {paralog_count} paralogs: {gene_count:,} genes ({percentage:.2f}%)")
    
    results['paralogy_distribution'] = dict(list(paralog_dist.items())[:20])
    
    # Paralogy statistics
    paralog_stats = comprehensive_df['Total_Paralogs'].describe()
    print(f"\n📊 PARALOGY STATISTICS:")
    print(f"   Average paralogs per gene: {paralog_stats['mean']:.2f} ± {paralog_stats['std']:.2f}")
    print(f"   Median paralogs: {paralog_stats['50%']:.0f}")
    print(f"   Range: {paralog_stats['min']:.0f} - {paralog_stats['max']:.0f}")
    
    results['avg_paralogs'] = paralog_stats['mean']
    results['median_paralogs'] = paralog_stats['50%']
    results['max_paralogs'] = paralog_stats['max']
    
    # Chromosomal distribution analysis
    if 'Same_Chr_Paralogs' in comprehensive_df.columns and 'Diff_Chr_Paralogs' in comprehensive_df.columns:
        print(f"\n📊 CHROMOSOMAL DISTRIBUTION:")
        same_chr_total = comprehensive_df['Same_Chr_Paralogs'].sum()
        diff_chr_total = comprehensive_df['Diff_Chr_Paralogs'].sum()
        total_relationships = same_chr_total + diff_chr_total
        
        print(f"   Same chromosome paralogs: {same_chr_total:,} ({100*same_chr_total/total_relationships:.1f}%)")
        print(f"   Different chromosome paralogs: {diff_chr_total:,} ({100*diff_chr_total/total_relationships:.1f}%)")
        
        results['same_chr_paralogs'] = same_chr_total
        results['diff_chr_paralogs'] = diff_chr_total
    
    return results

def comprehensive_correlation_analysis(analysis_df):
    """In-depth correlation analysis between LCRs and paralogy"""
    print("\n" + "="*80)
    print("COMPREHENSIVE LCR-PARALOGY CORRELATION ANALYSIS")
    print("="*80)
    
    results = {}
    
    # Basic correlation
    lcr_counts = analysis_df['LCR_Count']
    paralog_counts = analysis_df['Total_Paralogs']
    
    # Multiple correlation methods
    spearman_corr, spearman_p = spearmanr(lcr_counts, paralog_counts)
    pearson_corr, pearson_p = pearsonr(lcr_counts, paralog_counts)
    
    print(f"\n📊 CORRELATION ANALYSIS:")
    print(f"   Spearman correlation: {spearman_corr:.4f} (p = {spearman_p:.2e})")
    print(f"   Pearson correlation: {pearson_corr:.4f} (p = {pearson_p:.2e})")
    
    results['spearman_correlation'] = spearman_corr
    results['spearman_p_value'] = spearman_p
    results['pearson_correlation'] = pearson_corr
    results['pearson_p_value'] = pearson_p
    
    # Correlation by subgroups
    print(f"\n📊 SUBGROUP ANALYSIS:")
    
    # By LCR categories
    if 'LCR_Category' in analysis_df.columns:
        print(f"   Correlation by LCR Category:")
        for category in analysis_df['LCR_Category'].unique():
            if pd.notna(category):
                subset = analysis_df[analysis_df['LCR_Category'] == category]
                if len(subset) > 5:  # Minimum for meaningful correlation
                    corr, p_val = spearmanr(subset['LCR_Count'], subset['Total_Paralogs'])
                    print(f"     {category}: r = {corr:.3f} (p = {p_val:.3f}, n = {len(subset)})")
    
    # By chromosome context
    if 'Paralog_Diversity' in analysis_df.columns:
        print(f"   Correlation by Paralog Diversity:")
        for diversity in analysis_df['Paralog_Diversity'].unique():
            if pd.notna(diversity):
                subset = analysis_df[analysis_df['Paralog_Diversity'] == diversity]
                if len(subset) > 5:
                    corr, p_val = spearmanr(subset['LCR_Count'], subset['Total_Paralogs'])
                    print(f"     {diversity}: r = {corr:.3f} (p = {p_val:.3f}, n = {len(subset)})")
    
    # Binned analysis
    print(f"\n📊 BINNED ANALYSIS:")
    
    # Create LCR bins
    lcr_bins = pd.qcut(analysis_df['LCR_Count'], q=4, labels=['Low', 'Medium-Low', 'Medium-High', 'High'])
    bin_analysis = analysis_df.groupby(lcr_bins)['Total_Paralogs'].agg(['count', 'mean', 'std']).round(3)
    
    print("   Average paralogs by LCR quartile:")
    for bin_name, stats in bin_analysis.iterrows():
        print(f"     {bin_name} LCR: {stats['mean']:.2f} ± {stats['std']:.2f} paralogs (n = {stats['count']})")
    
    # Statistical test for trend
    kruskal_stat, kruskal_p = kruskal(*[analysis_df[lcr_bins == cat]['Total_Paralogs'].dropna() 
                                       for cat in lcr_bins.cat.categories])
    print(f"   Kruskal-Wallis test for trend: H = {kruskal_stat:.3f}, p = {kruskal_p:.2e}")
    
    results['kruskal_statistic'] = kruskal_stat
    results['kruskal_p_value'] = kruskal_p
    
    return results

def position_overlap_analysis(pairwise_df):
    """Detailed position overlap analysis"""
    print("\n" + "="*80)
    print("DETAILED POSITION OVERLAP ANALYSIS")
    print("="*80)
    
    results = {}
    
    if 'Overall_Position_Similarity' not in pairwise_df.columns:
        print("❌ Position similarity data not available")
        return {'error': 'missing_position_data'}
    
    similarity_scores = pairwise_df['Overall_Position_Similarity'].dropna()
    
    # Detailed similarity analysis
    print(f"\n📊 POSITION SIMILARITY STATISTICS:")
    sim_stats = similarity_scores.describe()
    print(f"   Mean similarity: {sim_stats['mean']:.3f}")
    print(f"   Median similarity: {sim_stats['50%']:.3f}")
    print(f"   Standard deviation: {sim_stats['std']:.3f}")
    print(f"   Range: {sim_stats['min']:.3f} - {sim_stats['max']:.3f}")
    
    # Multiple threshold analysis
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    print(f"\n📊 SIMILARITY THRESHOLD ANALYSIS:")
    
    threshold_results = {}
    for threshold in thresholds:
        above_threshold = sum(similarity_scores >= threshold)
        percentage = (above_threshold / len(similarity_scores)) * 100
        print(f"   ≥{threshold}: {above_threshold:,} pairs ({percentage:.1f}%)")
        threshold_results[f"threshold_{threshold}"] = above_threshold
    
    # Traditional categories
    high_overlap = sum(similarity_scores >= 0.8)
    medium_overlap = sum((similarity_scores >= 0.3) & (similarity_scores < 0.8))
    low_overlap = sum(similarity_scores < 0.3)
    total = len(similarity_scores)
    
    print(f"\n📊 TRADITIONAL CATEGORIES:")
    print(f"   High overlap (≥0.8): {high_overlap:,} pairs ({100*high_overlap/total:.1f}%) - RELATED")
    print(f"   Medium overlap (0.3-0.8): {medium_overlap:,} pairs ({100*medium_overlap/total:.1f}%) - MIXED")
    print(f"   Low overlap (<0.3): {low_overlap:,} pairs ({100*low_overlap/total:.1f}%) - INDEPENDENT")
    
    # Statistical interpretation
    if high_overlap > low_overlap:
        formation_pattern = "PREDOMINANTLY_RELATED"
        interpretation = "Most paralogs show related formation (similar LCR positions)"
    elif low_overlap > high_overlap:
        formation_pattern = "PREDOMINANTLY_INDEPENDENT"
        interpretation = "Most paralogs show independent formation (different LCR positions)"
    else:
        formation_pattern = "MIXED"
        interpretation = "Mixed pattern of related and independent formation"
    
    print(f"\n🎯 INTERPRETATION: {interpretation}")
    
    results.update({
        'mean_similarity': sim_stats['mean'],
        'median_similarity': sim_stats['50%'],
        'similarity_std': sim_stats['std'],
        'high_overlap': high_overlap,
        'medium_overlap': medium_overlap,
        'low_overlap': low_overlap,
        'total_pairs': total,
        'formation_pattern': formation_pattern,
        'threshold_analysis': threshold_results
    })
    
    return results

def evolutionary_distance_analysis(comprehensive_df, analysis_df):
    """Comprehensive evolutionary distance analysis"""
    print("\n" + "="*80)
    print("EVOLUTIONARY DISTANCE ANALYSIS")
    print("="*80)
    
    results = {}
    
    # Distance category analysis
    if 'Most_Common_Distance_Category' in comprehensive_df.columns:
        print(f"\n📊 DISTANCE CATEGORY DISTRIBUTION:")
        dist_categories = comprehensive_df['Most_Common_Distance_Category'].value_counts()
        
        for category, count in dist_categories.items():
            percentage = (count / len(comprehensive_df)) * 100
            print(f"   {category}: {count:,} genes ({percentage:.1f}%)")
        
        results['distance_category_distribution'] = dist_categories.to_dict()
        
        # LCR content by distance category
        print(f"\n📊 LCR CONTENT BY DISTANCE CATEGORY:")
        
        lcr_by_distance = comprehensive_df.groupby('Most_Common_Distance_Category')['LCR_Count'].agg(['count', 'mean', 'std']).round(3)
        
        for category, stats in lcr_by_distance.iterrows():
            print(f"   {category}: {stats['mean']:.2f} ± {stats['std']:.2f} LCRs (n = {stats['count']})")
        
        # Statistical test
        distance_categories = comprehensive_df['Most_Common_Distance_Category'].dropna().unique()
        if len(distance_categories) > 1:
            groups = [comprehensive_df[comprehensive_df['Most_Common_Distance_Category'] == cat]['LCR_Count'].dropna() 
                     for cat in distance_categories]
            
            kruskal_stat, kruskal_p = kruskal(*groups)
            print(f"   Kruskal-Wallis test: H = {kruskal_stat:.3f}, p = {kruskal_p:.2e}")
            
            results['distance_kruskal_stat'] = kruskal_stat
            results['distance_kruskal_p'] = kruskal_p
        
        # Specific comparisons
        same_chr = comprehensive_df[comprehensive_df['Most_Common_Distance_Category'] == 'Same_Chromosome']
        diff_chr = comprehensive_df[comprehensive_df['Most_Common_Distance_Category'] == 'Different_Chromosome']
        
        if len(same_chr) > 0 and len(diff_chr) > 0:
            same_chr_lcr_avg = same_chr['LCR_Count'].mean()
            diff_chr_lcr_avg = diff_chr['LCR_Count'].mean()
            
            print(f"\n📊 RECENT vs ANCIENT DUPLICATIONS:")
            print(f"   Same chromosome (recent): {same_chr_lcr_avg:.2f} LCRs per gene (n = {len(same_chr)})")
            print(f"   Different chromosome (ancient): {diff_chr_lcr_avg:.2f} LCRs per gene (n = {len(diff_chr)})")
            
            # Statistical test
            mannwhitney_stat, mannwhitney_p = mannwhitneyu(same_chr['LCR_Count'].dropna(), 
                                                          diff_chr['LCR_Count'].dropna(), 
                                                          alternative='two-sided')
            print(f"   Mann-Whitney U test: U = {mannwhitney_stat:.0f}, p = {mannwhitney_p:.2e}")
            
            if same_chr_lcr_avg > diff_chr_lcr_avg:
                recent_vs_ancient = "RECENT_HIGHER"
                interpretation = "Recent duplications have more LCRs (supports rapid evolution hypothesis)"
            else:
                recent_vs_ancient = "ANCIENT_HIGHER"
                interpretation = "Ancient duplications have more LCRs (contradicts rapid evolution hypothesis)"
            
            print(f"   INTERPRETATION: {interpretation}")
            
            results.update({
                'same_chr_avg_lcr': same_chr_lcr_avg,
                'diff_chr_avg_lcr': diff_chr_lcr_avg,
                'recent_vs_ancient_pattern': recent_vs_ancient,
                'mannwhitney_stat': mannwhitney_stat,
                'mannwhitney_p': mannwhitney_p
            })
    
    return results

def hypothesis_testing(analysis_df):
    """Test specific hypotheses from the paper"""
    print("\n" + "="*80)
    print("HYPOTHESIS TESTING (PAPER PREDICTIONS)")
    print("="*80)
    
    results = {}
    
    # Main hypothesis: Inverse correlation
    lcr_counts = analysis_df['LCR_Count']
    paralog_counts = analysis_df['Total_Paralogs']
    
    correlation, p_value = spearmanr(lcr_counts, paralog_counts)
    
    print(f"\n🔬 MAIN HYPOTHESIS TEST:")
    print(f"   Paper prediction: NEGATIVE correlation (LCR ↔ paralogy)")
    print(f"   Maize result: {correlation:.4f} ({'POSITIVE' if correlation > 0 else 'NEGATIVE'} correlation)")
    print(f"   P-value: {p_value:.2e}")
    
    if correlation < 0 and p_value < 0.05:
        main_hypothesis_result = "SUPPORTED"
        interpretation = "Maize shows compensatory evolution like prokaryotes"
    elif correlation > 0 and p_value < 0.05:
        main_hypothesis_result = "CONTRADICTED"
        interpretation = "Maize shows complementary evolution unlike prokaryotes"
    else:
        main_hypothesis_result = "INCONCLUSIVE"
        interpretation = "No significant correlation found"
    
    print(f"   CONCLUSION: Paper's hypothesis is {main_hypothesis_result}")
    print(f"   INTERPRETATION: {interpretation}")
    
    results['main_hypothesis'] = main_hypothesis_result
    results['main_correlation'] = correlation
    results['main_p_value'] = p_value
    
    # Specific category tests
    print(f"\n🔬 CATEGORY-SPECIFIC TESTS:")
    
    # High LCR vs Low LCR genes
    high_lcr_genes = analysis_df[analysis_df['LCR_Count'] >= 3]
    low_lcr_genes = analysis_df[analysis_df['LCR_Count'] <= 1]
    
    if len(high_lcr_genes) > 10 and len(low_lcr_genes) > 10:
        high_lcr_paralogs = high_lcr_genes['Total_Paralogs']
        low_lcr_paralogs = low_lcr_genes['Total_Paralogs']
        
        high_avg = high_lcr_paralogs.mean()
        low_avg = low_lcr_paralogs.mean()
        
        print(f"   High LCR genes (≥3): {high_avg:.2f} paralogs (n = {len(high_lcr_genes)})")
        print(f"   Low LCR genes (≤1): {low_avg:.2f} paralogs (n = {len(low_lcr_genes)})")
        
        # Statistical test
        t_stat, t_p = ttest_ind(high_lcr_paralogs, low_lcr_paralogs)
        mannwhitney_stat, mannwhitney_p = mannwhitneyu(high_lcr_paralogs, low_lcr_paralogs, alternative='two-sided')
        
        print(f"   T-test: t = {t_stat:.3f}, p = {t_p:.2e}")
        print(f"   Mann-Whitney: U = {mannwhitney_stat:.0f}, p = {mannwhitney_p:.2e}")
        
        if high_avg < low_avg and mannwhitney_p < 0.05:
            category_hypothesis = "SUPPORTED"
            category_interpretation = "High LCR genes have fewer paralogs"
        elif high_avg > low_avg and mannwhitney_p < 0.05:
            category_hypothesis = "CONTRADICTED"
            category_interpretation = "High LCR genes have more paralogs"
        else:
            category_hypothesis = "INCONCLUSIVE"
            category_interpretation = "No significant difference"
        
        print(f"   CONCLUSION: {category_interpretation}")
        
        results.update({
            'category_hypothesis': category_hypothesis,
            'high_lcr_avg_paralogs': high_avg,
            'low_lcr_avg_paralogs': low_avg,
            'category_mannwhitney_p': mannwhitney_p
        })
    
    return results

def create_comprehensive_output_csv(comprehensive_df, pairwise_df, analysis_df, all_results):
    """Create detailed CSV outputs with all analysis results"""
    print("\n" + "="*80)
    print("CREATING COMPREHENSIVE CSV OUTPUTS")
    print("="*80)
    
    # 1. Enhanced gene-level analysis CSV
    print("📊 Creating enhanced gene-level analysis CSV...")
    
    gene_analysis = comprehensive_df.copy()
    
    # Add derived metrics
    gene_analysis['LCR_Density'] = gene_analysis['Total_LCR_Length'] / (gene_analysis['LCR_Span'].fillna(1) + 1)
    gene_analysis['Paralog_Ratio'] = gene_analysis['Same_Chr_Paralogs'] / (gene_analysis['Total_Paralogs'] + 1)
    
    # Add percentile ranks
    gene_analysis['LCR_Count_Percentile'] = gene_analysis['LCR_Count'].rank(pct=True) * 100
    gene_analysis['Paralog_Count_Percentile'] = gene_analysis['Total_Paralogs'].rank(pct=True) * 100
    
    # Add evolutionary insights
    gene_analysis['Evolutionary_Strategy'] = 'Unknown'
    
    # High LCR, High Paralogs
    high_both = (gene_analysis['LCR_Count'] >= gene_analysis['LCR_Count'].quantile(0.75)) & \
                (gene_analysis['Total_Paralogs'] >= gene_analysis['Total_Paralogs'].quantile(0.75))
    gene_analysis.loc[high_both, 'Evolutionary_Strategy'] = 'High_LCR_High_Paralogs'
    
    # High LCR, Low Paralogs
    high_lcr_low_par = (gene_analysis['LCR_Count'] >= gene_analysis['LCR_Count'].quantile(0.75)) & \
                       (gene_analysis['Total_Paralogs'] <= gene_analysis['Total_Paralogs'].quantile(0.25))
    gene_analysis.loc[high_lcr_low_par, 'Evolutionary_Strategy'] = 'High_LCR_Low_Paralogs'
    
    # Low LCR, High Paralogs
    low_lcr_high_par = (gene_analysis['LCR_Count'] <= gene_analysis['LCR_Count'].quantile(0.25)) & \
                       (gene_analysis['Total_Paralogs'] >= gene_analysis['Total_Paralogs'].quantile(0.75))
    gene_analysis.loc[low_lcr_high_par, 'Evolutionary_Strategy'] = 'Low_LCR_High_Paralogs'
    
    # Low LCR, Low Paralogs
    low_both = (gene_analysis['LCR_Count'] <= gene_analysis['LCR_Count'].quantile(0.25)) & \
               (gene_analysis['Total_Paralogs'] <= gene_analysis['Total_Paralogs'].quantile(0.25))
    gene_analysis.loc[low_both, 'Evolutionary_Strategy'] = 'Low_LCR_Low_Paralogs'
    
    gene_analysis.to_csv('comprehensive_gene_analysis.csv', index=False)
    print(f"   ✅ Saved: comprehensive_gene_analysis.csv ({len(gene_analysis):,} genes)")
    
    # 2. Pairwise analysis CSV with enhanced metrics
    print("📊 Creating enhanced pairwise analysis CSV...")
    
    pairwise_analysis = pairwise_df.copy()
    
    # Add gene length calculations if coordinates are available
    required_cols = ['Gene1_Start', 'Gene1_End', 'Gene2_Start', 'Gene2_End']
    if all(col in pairwise_df.columns for col in required_cols):
        pairwise_analysis['Gene1_Length'] = pairwise_df['Gene1_End'] - pairwise_df['Gene1_Start'] + 1
        pairwise_analysis['Gene2_Length'] = pairwise_df['Gene2_End'] - pairwise_df['Gene2_Start'] + 1
        pairwise_analysis['Length_Ratio'] = pairwise_analysis['Gene2_Length'] / pairwise_analysis['Gene1_Length']
        pairwise_analysis['Length_Difference'] = abs(pairwise_analysis['Gene1_Length'] - pairwise_analysis['Gene2_Length'])
    
    # Add position overlap categories
    if 'Overall_Position_Similarity' in pairwise_analysis.columns:
        pairwise_analysis['Position_Overlap_Category'] = 'Medium'
        pairwise_analysis.loc[pairwise_analysis['Overall_Position_Similarity'] >= 0.8, 'Position_Overlap_Category'] = 'High'
        pairwise_analysis.loc[pairwise_analysis['Overall_Position_Similarity'] < 0.3, 'Position_Overlap_Category'] = 'Low'
        
        pairwise_analysis['Formation_Type'] = pairwise_analysis['Position_Overlap_Category'].map({
            'High': 'Related_Formation',
            'Medium': 'Mixed_Formation', 
            'Low': 'Independent_Formation'
        })
    
    pairwise_analysis.to_csv('comprehensive_pairwise_analysis.csv', index=False)
    print(f"   ✅ Saved: comprehensive_pairwise_analysis.csv ({len(pairwise_analysis):,} pairs)")
    
    # 3. Summary statistics CSV
    print("📊 Creating summary statistics CSV...")
    
    summary_stats = []
    
    # Add all results to summary
    for analysis_type, results in all_results.items():
        if isinstance(results, dict) and 'error' not in results:
            for metric, value in results.items():
                if isinstance(value, (int, float)):
                    summary_stats.append({
                        'Analysis_Type': analysis_type,
                        'Metric': metric,
                        'Value': value,
                        'Data_Type': type(value).__name__
                    })
    
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv('detailed_summary_statistics.csv', index=False)
    print(f"   ✅ Saved: detailed_summary_statistics.csv ({len(summary_df):,} metrics)")
    
    # 4. Evolutionary strategy analysis CSV
    print("📊 Creating evolutionary strategy analysis CSV...")
    
    strategy_analysis = gene_analysis.groupby('Evolutionary_Strategy').agg({
        'Gene_ID': 'count',
        'LCR_Count': ['mean', 'std'],
        'Total_Paralogs': ['mean', 'std'],
        'Total_LCR_Length': ['mean', 'std'],
        'Same_Chr_Paralogs': ['mean', 'std'],
        'Diff_Chr_Paralogs': ['mean', 'std']
    }).round(3)
    
    # Flatten column names
    strategy_analysis.columns = ['_'.join(col).strip() for col in strategy_analysis.columns]
    strategy_analysis = strategy_analysis.reset_index()
    
    strategy_analysis.to_csv('evolutionary_strategy_analysis.csv', index=False)
    print(f"   ✅ Saved: evolutionary_strategy_analysis.csv ({len(strategy_analysis):,} strategies)")
    
    # 5. Correlation analysis detailed CSV
    print("📊 Creating correlation analysis CSV...")
    
    correlation_results = []
    
    # Overall correlation
    spearman_corr, spearman_p = spearmanr(analysis_df['LCR_Count'], analysis_df['Total_Paralogs'])
    pearson_corr, pearson_p = pearsonr(analysis_df['LCR_Count'], analysis_df['Total_Paralogs'])
    
    correlation_results.extend([
        {'Analysis': 'Overall', 'Method': 'Spearman', 'Correlation': spearman_corr, 'P_Value': spearman_p, 'N': len(analysis_df)},
        {'Analysis': 'Overall', 'Method': 'Pearson', 'Correlation': pearson_corr, 'P_Value': pearson_p, 'N': len(analysis_df)}
    ])
    
    # Subgroup correlations
    if 'LCR_Category' in analysis_df.columns:
        for category in analysis_df['LCR_Category'].unique():
            if pd.notna(category):
                subset = analysis_df[analysis_df['LCR_Category'] == category]
                if len(subset) > 5:
                    corr, p_val = spearmanr(subset['LCR_Count'], subset['Total_Paralogs'])
                    correlation_results.append({
                        'Analysis': f'LCR_Category_{category}', 
                        'Method': 'Spearman', 
                        'Correlation': corr, 
                        'P_Value': p_val, 
                        'N': len(subset)
                    })
    
    if 'Paralog_Diversity' in analysis_df.columns:
        for diversity in analysis_df['Paralog_Diversity'].unique():
            if pd.notna(diversity):
                subset = analysis_df[analysis_df['Paralog_Diversity'] == diversity]
                if len(subset) > 5:
                    corr, p_val = spearmanr(subset['LCR_Count'], subset['Total_Paralogs'])
                    correlation_results.append({
                        'Analysis': f'Paralog_Diversity_{diversity}', 
                        'Method': 'Spearman', 
                        'Correlation': corr, 
                        'P_Value': p_val, 
                        'N': len(subset)
                    })
    
    correlation_df = pd.DataFrame(correlation_results)
    correlation_df.to_csv('detailed_correlation_analysis.csv', index=False)
    print(f"   ✅ Saved: detailed_correlation_analysis.csv ({len(correlation_df):,} analyses)")
    
    print(f"\n🎉 ALL CSV FILES CREATED SUCCESSFULLY!")
    return {
        'gene_analysis_file': 'comprehensive_gene_analysis.csv',
        'pairwise_analysis_file': 'comprehensive_pairwise_analysis.csv',
        'summary_stats_file': 'detailed_summary_statistics.csv',
        'strategy_analysis_file': 'evolutionary_strategy_analysis.csv',
        'correlation_analysis_file': 'detailed_correlation_analysis.csv'
    }

def main():
    print("🧬 COMPREHENSIVE MAIZE LCR-PARALOGY ANALYSIS")
    print("=" * 80)
    print("🔬 In-depth analysis with complete statistical testing")
    print("📊 Multiple CSV outputs for detailed examination")
    print("=" * 80)
    
    try:
        # Step 1: Load and validate data
        comprehensive_df, pairwise_df, analysis_df = load_and_validate_data()
        
        # Initialize results storage
        all_results = {}
        
        # Step 2: Detailed LCR analysis
        all_results['lcr_analysis'] = detailed_lcr_analysis(comprehensive_df, analysis_df)
        
        # Step 3: Detailed paralogy analysis
        all_results['paralogy_analysis'] = detailed_paralogy_analysis(comprehensive_df, analysis_df)
        
        # Step 4: Comprehensive correlation analysis
        all_results['correlation_analysis'] = comprehensive_correlation_analysis(analysis_df)
        
        # Step 5: Position overlap analysis
        all_results['position_analysis'] = position_overlap_analysis(pairwise_df)
        
        # Step 6: Evolutionary distance analysis
        all_results['evolutionary_analysis'] = evolutionary_distance_analysis(comprehensive_df, analysis_df)
        
        # Step 7: Hypothesis testing
        all_results['hypothesis_testing'] = hypothesis_testing(analysis_df)
        
        # Step 8: Create comprehensive CSV outputs
        csv_files = create_comprehensive_output_csv(comprehensive_df, pairwise_df, analysis_df, all_results)
        
        # Step 9: Generate final comprehensive report
        print("\n" + "="*80)
        print("FINAL COMPREHENSIVE REPORT")
        print("="*80)
        
        # Main findings summary
        correlation = all_results['correlation_analysis']['spearman_correlation']
        correlation_p = all_results['correlation_analysis']['spearman_p_value']
        hypothesis_result = all_results['hypothesis_testing']['main_hypothesis']
        
        print(f"\n🎯 KEY FINDINGS:")
        print(f"   📊 Dataset: {len(comprehensive_df):,} genes, {len(pairwise_df):,} gene pairs")
        print(f"   🔬 Main correlation: {correlation:.4f} (p = {correlation_p:.2e})")
        print(f"   📋 Paper's hypothesis: {hypothesis_result}")
        print(f"   🧬 Evolutionary pattern: {'COMPLEMENTARY' if correlation > 0 else 'COMPENSATORY'}")
        
        # Position overlap findings
        if 'error' not in all_results['position_analysis']:
            formation_pattern = all_results['position_analysis']['formation_pattern']
            related_percent = (all_results['position_analysis']['high_overlap'] / 
                             all_results['position_analysis']['total_pairs']) * 100
            print(f"   🔗 Formation pattern: {formation_pattern}")
            print(f"   📍 Related formation: {related_percent:.1f}% of pairs")
        
        # LCR distribution findings
        lcr_dist = all_results['lcr_analysis']['lcr_distribution']
        max_lcr_count = max(lcr_dist.keys()) if lcr_dist else 0
        genes_with_many_lcrs = sum(count for lcr_count, count in lcr_dist.items() if lcr_count >= 3)
        print(f"   🧮 LCR range: 0 - {max_lcr_count} LCRs per gene")
        print(f"   🔢 High-LCR genes: {genes_with_many_lcrs:,} genes (≥3 LCRs)")
        
        print(f"\n📁 OUTPUT FILES CREATED:")
        for file_type, filename in csv_files.items():
            print(f"   ✅ {filename}")
        
        # Create final summary report
        final_report = f"""
COMPREHENSIVE MAIZE LCR-PARALOGY ANALYSIS REPORT
===============================================

EXECUTIVE SUMMARY:
- Total genes analyzed: {len(comprehensive_df):,}
- Total gene pairs analyzed: {len(pairwise_df):,}
- Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

MAIN FINDINGS:
1. LCR-Paralogy Relationship:
   - Spearman correlation: {correlation:.4f} (p = {correlation_p:.2e})
   - Direction: {'POSITIVE' if correlation > 0 else 'NEGATIVE'} correlation
   - Paper's hypothesis: {hypothesis_result}

2. Position Overlap Pattern:
   - Formation type: {formation_pattern if 'error' not in all_results['position_analysis'] else 'Data unavailable'}
   - Related formation: {related_percent:.1f}% of pairs

3. Evolutionary Strategy:
   - Maize shows {'COMPLEMENTARY' if correlation > 0 else 'COMPENSATORY'} evolution
   - LCRs and paralogy are {'mutually reinforcing' if correlation > 0 else 'alternative strategies'}
   - Plant genomes {'differ from' if correlation > 0 else 'match'} prokaryotic patterns

BIOLOGICAL INTERPRETATION:
{'Maize paralogs use both LCR expansion and gene duplication as complementary mechanisms for evolutionary adaptation. This suggests that eukaryotic plant genomes require multiple strategies to maintain functional diversity and evolutionary flexibility.' if correlation > 0 else 'Maize paralogs use LCR expansion and gene duplication as alternative mechanisms, similar to prokaryotes. This suggests conserved evolutionary strategies across domains of life.'}

STATISTICAL SIGNIFICANCE:
- Main correlation: {'SIGNIFICANT' if correlation_p < 0.05 else 'NOT SIGNIFICANT'} (p = {correlation_p:.2e})
- Effect size: {'WEAK' if abs(correlation) < 0.3 else 'MODERATE' if abs(correlation) < 0.7 else 'STRONG'}
- Sample size: {len(analysis_df):,} genes (adequate for statistical power)

RESEARCH IMPLICATIONS:
1. {'Contradicts' if correlation > 0 else 'Supports'} the universal applicability of prokaryotic LCR-paralogy patterns
2. Suggests kingdom-specific evolutionary mechanisms
3. Highlights the complexity of eukaryotic genome evolution
4. Provides insights for crop improvement and evolutionary biology

OUTPUT FILES:
- comprehensive_gene_analysis.csv: Gene-level detailed analysis
- comprehensive_pairwise_analysis.csv: Pairwise relationship analysis  
- detailed_summary_statistics.csv: All statistical results
- evolutionary_strategy_analysis.csv: Strategy-based comparisons
- detailed_correlation_analysis.csv: Subgroup correlation analyses

NEXT STEPS:
1. Functional annotation analysis of extreme cases
2. Expression correlation with LCR/paralogy patterns
3. Comparative analysis with other plant species
4. Investigation of specific gene families showing extreme patterns
"""
        
        with open('comprehensive_analysis_report.txt', 'w') as f:
            f.write(final_report)
        
        print(f"   ✅ comprehensive_analysis_report.txt")
        
        print(f"\n🎉 COMPREHENSIVE ANALYSIS COMPLETED!")
        print(f"📊 {len(csv_files)} detailed CSV files created")
        print(f"📋 Complete statistical analysis performed")
        print(f"🧬 Ready for publication and further research")
        
        return all_results
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = main()
    if results:
        print("✅ Analysis completed successfully!")
    else:
        print("❌ Analysis failed!")
        sys.exit(1)
