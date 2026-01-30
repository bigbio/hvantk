import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(14, 10))
ax.set_xlim(0, 14)
ax.set_ylim(0, 10)
ax.axis('off')

# Colors
color_vds = '#4A90A4'       # Blue-gray for VDS/input
color_hvantk = '#E8F4F8'    # Light blue for hvantk box
color_genomic = '#2E7D32'   # Green for genomic
color_transcr = '#1565C0'   # Blue for transcriptomic
color_proteomic = '#7B1FA2' # Purple for proteomic
color_psroc = '#F57C00'     # Orange for PSROC
color_enrichex = '#C62828'  # Red for EnrichEx
color_qc = '#455A64'        # Gray for QC
color_output = '#37474F'    # Dark gray for outputs

def draw_box(ax, x, y, w, h, color, text, fontsize=10, text_color='white', alpha=1.0, bold=False):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.15",
                          facecolor=color, edgecolor='none', alpha=alpha)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fontsize,
            color=text_color, weight=weight, wrap=True)

def draw_arrow(ax, start, end, color='#666666', style='->'):
    arrow = FancyArrowPatch(start, end, arrowstyle=style, color=color,
                            mutation_scale=15, lw=2, connectionstyle="arc3,rad=0")
    ax.add_patch(arrow)

# Title
ax.text(7, 9.6, 'hvantk Workflow', ha='center', va='center', fontsize=18, weight='bold', color='#1a1a1a')

# Input section - VDS
draw_box(ax, 0.5, 7.5, 2.5, 1.2, color_vds, 'VDS\n(from VDS Combiner)', fontsize=11, bold=True)
ax.text(1.75, 6.9, 'King et al., 2024', ha='center', va='center', fontsize=8, style='italic', color='#666')

# Main hvantk box (background)
hvantk_box = FancyBboxPatch((3.5, 1.0), 10, 7.5, boxstyle="round,pad=0.02,rounding_size=0.3",
                             facecolor=color_hvantk, edgecolor='#4A90A4', linewidth=3, alpha=0.3)
ax.add_patch(hvantk_box)
ax.text(8.5, 8.2, 'hvantk', ha='center', va='center', fontsize=16, weight='bold', color='#4A90A4')

# Arrow from VDS to hvantk
draw_arrow(ax, (3.0, 8.1), (3.8, 8.1))

# ============ Multi-omics Annotation Section ============
ax.text(6.0, 7.3, 'Multi-omics Annotation Framework', ha='center', va='center', 
        fontsize=12, weight='bold', color='#333')

# Genomic annotations
draw_box(ax, 4.0, 5.8, 1.8, 1.2, color_genomic, 'Genomic\n\nClinVar\ndbNSFP\ngnomAD', fontsize=8)

# Transcriptomic annotations  
draw_box(ax, 6.0, 5.8, 1.8, 1.2, color_transcr, 'Transcriptomic\n\nGTEx\nExpr. Atlas\nUCSC Cell', fontsize=8)

# Proteomic annotations
draw_box(ax, 8.0, 5.8, 1.8, 1.2, color_proteomic, 'Proteomic\n\nCPTAC\nINSIDER', fontsize=8)

# ============ Analysis Modules Section ============
ax.text(8.5, 4.8, 'Analysis Modules', ha='center', va='center', 
        fontsize=12, weight='bold', color='#333')

# PSROC module
draw_box(ax, 4.0, 3.2, 2.5, 1.2, color_psroc, 'PSROC\n\nROC Analysis\nThreshold Optimization', fontsize=9)

# EnrichEx module
draw_box(ax, 7.0, 3.2, 2.5, 1.2, color_enrichex, 'EnrichEx\n\nOverlap Enrichment\nBurden Testing', fontsize=9)

# QC module
draw_box(ax, 10.0, 3.2, 2.5, 1.2, color_qc, 'QC & Reporting\n\nSample/Variant QC\nHTML Reports', fontsize=9)

# ============ Outputs Section ============
ax.text(8.5, 2.2, 'Outputs', ha='center', va='center', 
        fontsize=12, weight='bold', color='#333')

# Output boxes
draw_box(ax, 4.0, 1.2, 2.2, 0.8, color_output, 'Annotated\nHail Tables', fontsize=9)
draw_box(ax, 6.5, 1.2, 2.2, 0.8, color_output, 'ROC Curves\n& Metrics', fontsize=9)
draw_box(ax, 9.0, 1.2, 2.2, 0.8, color_output, 'Enrichment\nResults', fontsize=9)
draw_box(ax, 11.5, 1.2, 1.8, 0.8, color_output, 'QC\nReports', fontsize=9)

# ============ Arrows ============
# From VDS input into annotation framework (implied by box)
# Vertical arrows from annotation to analysis
draw_arrow(ax, (5.9, 5.8), (5.25, 4.4), color='#888')
draw_arrow(ax, (6.9, 5.8), (8.25, 4.4), color='#888')
draw_arrow(ax, (8.9, 5.8), (11.25, 4.4), color='#888')

# Arrows from modules to outputs
draw_arrow(ax, (5.1, 5.8), (5.1, 2.0), color='#888')
draw_arrow(ax, (5.25, 3.2), (5.1, 2.0), color='#888')
draw_arrow(ax, (8.25, 3.2), (7.6, 2.0), color='#888')
draw_arrow(ax, (8.25, 3.2), (10.1, 2.0), color='#888')
draw_arrow(ax, (11.25, 3.2), (12.2, 2.0), color='#888')

# Legend
legend_y = 0.3
ax.add_patch(FancyBboxPatch((0.5, legend_y), 0.3, 0.3, boxstyle="round,pad=0.01", facecolor=color_genomic, edgecolor='none'))
ax.text(1.0, legend_y + 0.15, 'Genomic', fontsize=8, va='center')

ax.add_patch(FancyBboxPatch((2.0, legend_y), 0.3, 0.3, boxstyle="round,pad=0.01", facecolor=color_transcr, edgecolor='none'))
ax.text(2.5, legend_y + 0.15, 'Transcriptomic', fontsize=8, va='center')

ax.add_patch(FancyBboxPatch((4.0, legend_y), 0.3, 0.3, boxstyle="round,pad=0.01", facecolor=color_proteomic, edgecolor='none'))
ax.text(4.5, legend_y + 0.15, 'Proteomic', fontsize=8, va='center')

ax.add_patch(FancyBboxPatch((5.8, legend_y), 0.3, 0.3, boxstyle="round,pad=0.01", facecolor=color_psroc, edgecolor='none'))
ax.text(6.3, legend_y + 0.15, 'PSROC', fontsize=8, va='center')

ax.add_patch(FancyBboxPatch((7.5, legend_y), 0.3, 0.3, boxstyle="round,pad=0.01", facecolor=color_enrichex, edgecolor='none'))
ax.text(8.0, legend_y + 0.15, 'EnrichEx', fontsize=8, va='center')

ax.add_patch(FancyBboxPatch((9.3, legend_y), 0.3, 0.3, boxstyle="round,pad=0.01", facecolor=color_qc, edgecolor='none'))
ax.text(9.8, legend_y + 0.15, 'QC', fontsize=8, va='center')

plt.tight_layout()
plt.savefig('/Users/yperez/work/hvantk/paper/figures/fig1_workflow.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.savefig('/Users/yperez/work/hvantk/paper/figures/fig1_workflow.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("Workflow figure saved!")
