#!/usr/bin/env python3
"""
Web Registry Generator for HVANTK Dataset Validation Registry

This module generates static web interfaces for the dataset validation registry,
including interactive dashboards, charts, and detailed status pages.
"""
import json
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any
from collections import defaultdict, Counter

from .config import RegistryConfig


class WebRegistryGenerator:
    """Generator for web-based dataset validation registry interface."""

    def __init__(self, config: RegistryConfig = None):
        """Initialize the web generator with configuration."""
        self.config = config or RegistryConfig()

    def generate_status_badge(self, status: str) -> str:
        """Generate HTML badge for dataset status."""
        color = self.config.chart_colors.get(status, '#6c757d')
        return f'<span class="badge" style="background-color: {color}; color: white; padding: 4px 8px; border-radius: 4px; font-size: 0.8em;">{status.replace("_", " ").title()}</span>'

    def generate_dataset_table(self, datasets: Dict[str, Any], dataset_type: str) -> str:
        """Generate HTML table for datasets."""
        if not datasets:
            return f"<p>No {dataset_type} datasets found.</p>"

        html = f"""
        <table class="dataset-table">
            <thead>
                <tr>
                    <th>Dataset ID</th>
                    <th>Status</th>
                    <th>Last Updated</th>
                    <th>Error Message</th>
                    <th>Details</th>
                </tr>
            </thead>
            <tbody>
        """

        # Sort datasets by status (passed first) then by ID
        sorted_datasets = sorted(
            datasets.items(),
            key=lambda x: (x[1]['status'] not in ['tier3_passed', 'tier2_passed'], x[0])
        )

        for dataset_id, data in sorted_datasets:
            status_badge = self.generate_status_badge(data['status'])
            timestamp = datetime.fromisoformat(data['timestamp']).strftime('%Y-%m-%d %H:%M')

            # Fix error message handling for None values
            error_msg = data.get('error_message') or ''
            if len(error_msg) > 100:
                error_msg = error_msg[:100] + '...'

            html += f"""
                <tr class="{'success' if data['status'] in ['tier3_passed', 'tier2_passed'] else 'failed'}">
                    <td><strong>{dataset_id}</strong></td>
                    <td>{status_badge}</td>
                    <td>{timestamp}</td>
                    <td class="error-message">{error_msg}</td>
                    <td><a href="#" onclick="showDetails('{dataset_id}')">View Details</a></td>
                </tr>
            """

        html += """
            </tbody>
        </table>
        """
        return html

    def generate_summary_stats(self, registry: Dict[str, Any]) -> Dict[str, Any]:
        """Generate summary statistics."""
        total_datasets = len(registry)
        status_counts = Counter(data['status'] for data in registry.values())
        type_counts = Counter(data['dataset_type'] for data in registry.values())

        successful = sum(count for status, count in status_counts.items()
                        if status in ['tier3_passed', 'tier2_passed'])

        return {
            'total': total_datasets,
            'successful': successful,
            'failed': total_datasets - successful,
            'success_rate': round((successful / total_datasets * 100) if total_datasets > 0 else 0, 1),
            'status_counts': dict(status_counts),
            'type_counts': dict(type_counts),
            'last_updated': max((data['timestamp'] for data in registry.values()), default='Never')
        }

    def generate_index_html(self, registry: Dict[str, Any], output_dir: Path) -> None:
        """Generate the main index.html file."""
        stats = self.generate_summary_stats(registry)

        # Separate datasets by type
        ucsc_datasets = {k: v for k, v in registry.items() if v['dataset_type'] == 'ucsc'}
        atlas_datasets = {k: v for k, v in registry.items() if v['dataset_type'] == 'expression_atlas'}

        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.config.site_title}</title>
    <link rel="stylesheet" href="styles.css">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
</head>
<body>
    <header>
        <div class="container">
            <h1>🧬 {self.config.site_title}</h1>
            <p>{self.config.site_description}</p>
        </div>
    </header>

    <main class="container">
        <!-- Summary Statistics -->
        <section class="stats-section">
            <h2>📊 Summary Statistics</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>{stats['total']}</h3>
                    <p>Total Datasets</p>
                </div>
                <div class="stat-card success">
                    <h3>{stats['successful']}</h3>
                    <p>Validated Successfully</p>
                </div>
                <div class="stat-card failed">
                    <h3>{stats['failed']}</h3>
                    <p>Failed Validation</p>
                </div>
                <div class="stat-card">
                    <h3>{stats['success_rate']}%</h3>
                    <p>Success Rate</p>
                </div>
            </div>
            <p class="last-updated">Last updated: {datetime.fromisoformat(stats['last_updated']).strftime('%Y-%m-%d %H:%M UTC') if stats['last_updated'] != 'Never' else 'Never'}</p>
        </section>

        <!-- Validation Status Chart -->
        <section class="chart-section">
            <h2>📈 Validation Status Distribution</h2>
            <canvas id="statusChart" width="400" height="200"></canvas>
        </section>

        <!-- UCSC Datasets -->
        <section class="datasets-section">
            <h2>🔬 UCSC Cell Browser Datasets ({len(ucsc_datasets)} datasets)</h2>
            {self.generate_dataset_table(ucsc_datasets, 'UCSC')}
        </section>

        <!-- Expression Atlas Datasets -->
        <section class="datasets-section">
            <h2>🧪 Expression Atlas Datasets ({len(atlas_datasets)} datasets)</h2>
            {self.generate_dataset_table(atlas_datasets, 'Expression Atlas')}
        </section>

        <!-- Legend -->
        <section class="legend-section">
            <h2>📖 Status Legend</h2>
            <div class="legend">
                <div class="legend-item">
                    {self.generate_status_badge('tier3_passed')}
                    <span>Tier 3 Passed - Full dataset validated successfully</span>
                </div>
                <div class="legend-item">
                    {self.generate_status_badge('tier2_passed')}
                    <span>Tier 2 Passed - Sample validation successful</span>
                </div>
                <div class="legend-item">
                    {self.generate_status_badge('tier1_passed')}
                    <span>Tier 1 Passed - Header validation successful</span>
                </div>
                <div class="legend-item">
                    {self.generate_status_badge('tier1_failed')}
                    <span>Failed - Header validation failed</span>
                </div>
                <div class="legend-item">
                    {self.generate_status_badge('download_failed')}
                    <span>Download Failed - Unable to download dataset</span>
                </div>
            </div>
        </section>
    </main>

    <!-- Dataset Details Modal -->
    <div id="detailsModal" class="modal">
        <div class="modal-content">
            <span class="close">&times;</span>
            <div id="modalBody"></div>
        </div>
    </div>

    <footer>
        <div class="container">
            <p>Generated automatically by <a href="https://github.com/your-org/pyvatk">HVANTK</a> | 
               <a href="https://github.com/your-org/pyvatk/actions">View Validation Runs</a></p>
        </div>
    </footer>

    <script>
        // Dataset registry data for JavaScript
        const registryData = {json.dumps(registry, indent=2)};
        
        // Chart.js configuration  
        const ctx = document.getElementById('statusChart').getContext('2d');
        const statusChart = new Chart(ctx, {{
            type: 'doughnut',
            data: {{
                labels: {list(stats['status_counts'].keys())},
                datasets: [{{
                    data: {list(stats['status_counts'].values())},
                    backgroundColor: {[self.config.chart_colors.get(status, '#6c757d') for status in stats['status_counts'].keys()]}
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{
                    legend: {{
                        position: 'bottom'
                    }}
                }}
            }}
        }});
        
        // Modal functionality
        const modal = document.getElementById('detailsModal');
        const span = document.getElementsByClassName('close')[0];
        
        function showDetails(datasetId) {{
            const dataset = registryData[datasetId];
            document.getElementById('modalBody').innerHTML = generateDetailView(datasetId, dataset);
            modal.style.display = 'block';
        }}
        
        span.onclick = function() {{
            modal.style.display = 'none';
        }}
        
        window.onclick = function(event) {{
            if (event.target == modal) {{
                modal.style.display = 'none';
            }}
        }}
        
        function generateDetailView(datasetId, dataset) {{
            return `
                <h2>${{datasetId}}</h2>
                <p><strong>Type:</strong> ${{dataset.dataset_type}}</p>
                <p><strong>Status:</strong> ${{dataset.status}}</p>
                <p><strong>Last Updated:</strong> ${{new Date(dataset.timestamp).toLocaleString()}}</p>
                ${{dataset.error_message ? `<p><strong>Error:</strong> ${{dataset.error_message}}</p>` : ''}}
                ${{dataset.tier1_details ? `<h3>Tier 1 Details</h3><pre>${{JSON.stringify(dataset.tier1_details, null, 2)}}</pre>` : ''}}
                ${{dataset.tier2_details ? `<h3>Tier 2 Details</h3><pre>${{JSON.stringify(dataset.tier2_details, null, 2)}}</pre>` : ''}}
                ${{dataset.tier3_details ? `<h3>Tier 3 Details</h3><pre>${{JSON.stringify(dataset.tier3_details, null, 2)}}</pre>` : ''}}
            `;
        }}
    </script>
</body>
</html>
        """

        with open(output_dir / 'index.html', 'w') as f:
            f.write(html)

    def generate_css(self, output_dir: Path) -> None:
        """Generate the CSS file."""
        css = f"""
/* HVANTK Dataset Registry Styles */
* {{
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}}

body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    line-height: 1.6;
    color: #333;
    background-color: #f8f9fa;
}}

.container {{
    max-width: 1200px;
    margin: 0 auto;
    padding: 0 20px;
}}

header {{
    background: linear-gradient(135deg, {self.config.theme_color} 0%, #764ba2 100%);
    color: white;
    padding: 2rem 0;
    text-align: center;
}}

header h1 {{
    font-size: 2.5rem;
    margin-bottom: 0.5rem;
}}

main {{
    padding: 2rem 0;
}}

.stats-section {{
    margin-bottom: 2rem;
}}

.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    margin: 1rem 0;
}}

.stat-card {{
    background: white;
    padding: 1.5rem;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    text-align: center;
    border-left: 4px solid {self.config.theme_color};
}}

.stat-card.success {{
    border-left-color: #28a745;
}}

.stat-card.failed {{
    border-left-color: #dc3545;
}}

.stat-card h3 {{
    font-size: 2rem;
    color: #333;
    margin-bottom: 0.5rem;
}}

.chart-section {{
    background: white;
    padding: 2rem;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 2rem;
}}

.chart-section canvas {{
    max-height: 400px;
}}

.datasets-section {{
    background: white;
    padding: 2rem;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 2rem;
}}

.dataset-table {{
    width: 100%;
    border-collapse: collapse;
    margin-top: 1rem;
}}

.dataset-table th,
.dataset-table td {{
    padding: 12px;
    text-align: left;
    border-bottom: 1px solid #ddd;
}}

.dataset-table th {{
    background-color: #f8f9fa;
    font-weight: 600;
}}

.dataset-table tr.success {{
    background-color: #f8fff9;
}}

.dataset-table tr.failed {{
    background-color: #fff8f8;
}}

.dataset-table tr:hover {{
    background-color: #f5f5f5;
}}

.error-message {{
    font-family: monospace;
    font-size: 0.9em;
    color: #666;
}}

.badge {{
    display: inline-block;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 0.8em;
    font-weight: 500;
}}

.legend-section {{
    background: white;
    padding: 2rem;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}}

.legend {{
    display: flex;
    flex-direction: column;
    gap: 1rem;
}}

.legend-item {{
    display: flex;
    align-items: center;
    gap: 1rem;
}}

.last-updated {{
    text-align: center;
    color: #666;
    font-style: italic;
    margin-top: 1rem;
}}

/* Modal Styles */
.modal {{
    display: none;
    position: fixed;
    z-index: 1000;
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    background-color: rgba(0,0,0,0.5);
}}

.modal-content {{
    background-color: white;
    margin: 5% auto;
    padding: 2rem;
    border-radius: 8px;
    width: 80%;
    max-width: 800px;
    max-height: 80%;
    overflow-y: auto;
}}

.close {{
    color: #aaa;
    float: right;
    font-size: 28px;
    font-weight: bold;
    cursor: pointer;
}}

.close:hover {{
    color: black;
}}

footer {{
    background: #333;
    color: white;
    text-align: center;
    padding: 1rem 0;
    margin-top: 2rem;
}}

footer a {{
    color: {self.config.theme_color};
    text-decoration: none;
}}

footer a:hover {{
    text-decoration: underline;
}}

/* Responsive Design */
@media (max-width: 768px) {{
    .container {{
        padding: 0 10px;
    }}
    
    header h1 {{
        font-size: 2rem;
    }}
    
    .stats-grid {{
        grid-template-columns: 1fr;
    }}
    
    .dataset-table {{
        font-size: 0.9em;
    }}
    
    .modal-content {{
        width: 95%;
        margin: 10% auto;
    }}
}}
"""

        with open(output_dir / 'styles.css', 'w') as f:
            f.write(css)

    def generate_web_registry(self, registry_file: str, output_dir: str) -> None:
        """Generate complete web registry from validation data."""
        # Load registry data
        with open(registry_file, 'r') as f:
            registry = json.load(f)

        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Generate web files
        self.generate_index_html(registry, output_path)
        self.generate_css(output_path)

        print(f"✅ Web registry generated successfully in {output_path}")
        return output_path
