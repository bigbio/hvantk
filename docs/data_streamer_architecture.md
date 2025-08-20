# Data Streamer Architecture Documentation

## Overview

The Data Streamer architecture provides a scalable, memory-efficient way to process large genomic datasets in chunks. This implementation focuses on the Clinvar training set generation use case but is designed to be extensible for other data processing pipelines.

## Architecture Components

### 1. Base Classes

- **`DataStreamer`**: Abstract base class defining the streaming interface
- **`HailDataStreamer`**: Specialized base for Hail-based genomic data processing
- **`StreamProcessor`**: Orchestrates multiple streamers in a pipeline

### 2. Clinvar Implementation

- **`ClinvarDataStreamer`**: Processes Clinvar VCF data in chunks
- **`ClinvarTrainingSetProcessor`**: Complete pipeline for training set generation

## Key Features

- **Memory Efficiency**: Processes data in configurable chunks (default: 10,000 variants)
- **Modular Design**: Each processing step is a separate, reusable streamer
- **Hail Integration**: Automatic Hail initialization and cleanup
- **Error Handling**: Proper setup/teardown with resource management
- **Extensibility**: Easy to add new streamers for different data sources

## Usage Examples

### Basic Usage

```python
from hvantk.utils.clinvar_streamer import create_clinvar_training_set_streamer

# Create processor
processor = create_clinvar_training_set_streamer(
    clinvar_path="./data/clinvar/clinvar_20220403.vcf.gz",
    output_dir="./data/training_set"
)

# Process data
training_set = processor.process()

if training_set:
    print(f"Generated {training_set.count()} training examples")
```

### Custom Configuration

```python
from hvantk.utils.clinvar_streamer import ClinvarDataStreamer, ClinvarTrainingSetProcessor

# Define custom CHD genes
chd_genes = {"GATA4", "NKX2-5", "TBX5", "NOTCH1", "CHD7"}

# Create custom processor
processor = ClinvarTrainingSetProcessor(
    clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
    chd_genes=chd_genes,
    output_dir="./output"
)

result = processor.process()
```

### Using Individual Streamers

```python
from hvantk.utils.clinvar_streamer import ClinvarDataStreamer

# Create streamer with custom chunk size
streamer = ClinvarDataStreamer(
    clinvar_path="./data/clinvar/clinvar.vcf.gz",
    chd_genes={"GATA4", "NKX2-5"},
    chunk_size=5000  # Smaller chunks
)

# Manual processing
streamer.setup()
try:
    for chunk in streamer.stream():
        # Process each chunk individually
        print(f"Processed chunk with {chunk.count()} variants")
finally:
    streamer.teardown()
```

## Configuration

### Clinvar Labels

The streamer uses predefined label sets for classification:

- **Pathogenic**: "Pathogenic/Likely_pathogenic", "Likely_pathogenic", "Pathogenic"
- **Benign**: "Benign/Likely_benign", "Likely_benign", "Benign"
- **CHD-specific**: "Congenital_heart_disease", "Congenital_heart_defect"

### Training Set Logic

1. **True Positives (TP)**: Variants that are:
   - Pathogenic AND in CHD genes, OR
   - Associated with CHD diseases
   
2. **True Negatives (TN)**: Variants that are:
   - Benign (regardless of gene)

3. **Filtering**: Only variants that are exclusively TP or TN (not both)

## Output

The processor generates:
- **Hail Table**: `ts.clinvar.ht` (native Hail format)
- **TSV file**: `ts.clinvar.ht.tsv` (human-readable export)

Output columns:
- `gene`: Gene symbol
- `rf_label`: Training label ("TP" or "TN")

## Performance Considerations

### Memory Usage
- Chunk size controls memory usage vs. processing efficiency
- Default 10,000 variants per chunk balances both
- Adjust based on available memory and dataset size

### Disk Usage
- Intermediate results are checkpointed to disk
- Final output includes both native Hail format and TSV export

## Extending the Architecture

### Adding New Streamers

```python
from hvantk.data.data_streamer import HailDataStreamer
import hail as hl

class MyCustomStreamer(HailDataStreamer):
    def __init__(self, data_path: str):
        super().__init__("MyCustomStreamer")
        self.data_path = data_path
        
    def setup(self):
        super().setup()
        # Load your data
        self.data = hl.import_table(self.data_path)
        
    def stream(self):
        # Implement chunking logic
        total = self.data.count()
        for i in range(0, total, self.chunk_size):
            chunk = self.data.head(self.chunk_size)
            yield self.process_chunk(chunk)
            
    def process_chunk(self, chunk):
        # Implement chunk processing
        return chunk.annotate(processed=True)
```

### Creating Multi-Step Pipelines

```python
from hvantk.data.data_streamer import StreamProcessor

# Create pipeline with multiple streamers
pipeline = StreamProcessor("MyPipeline")
pipeline.add_streamer(ClinvarDataStreamer(...))
pipeline.add_streamer(MyCustomStreamer(...))

# Process through entire pipeline
result = pipeline.process()
```

## Error Handling

The architecture includes robust error handling:
- Automatic Hail initialization/cleanup
- Proper resource management with try/finally blocks
- Detailed logging at each step
- Graceful handling of empty datasets

## Logging

Each streamer maintains its own logger:
```python
import logging

# Configure logging level
logging.basicConfig(level=logging.INFO)

# Streamers will log progress and statistics
processor = create_clinvar_training_set_streamer(...)
result = processor.process()  # Logs progress automatically
```

## Migration from Original Script

The new streamer architecture maintains the same functionality as the original `generate_training_set.py` script but with these advantages:

1. **Better Memory Management**: Processes data in chunks
2. **Improved Modularity**: Each step is a separate, testable component
3. **Enhanced Reusability**: Streamers can be used in different combinations
4. **Better Error Handling**: Robust resource management
5. **Extensibility**: Easy to add new data sources and processing steps

To migrate existing code:
```python
# Old way
from generate_training_set import main
main()

# New way
from hvantk.utils.generate_training_set_streamed import main
main()
```

## Testing

Basic tests are included in `hvantk/tests/test_data_streamer.py`. Run with:
```bash
python -m hvantk.tests.test_data_streamer
```

## Next Steps

1. **Implement Real CHD Gene Loading**: Replace the placeholder `load_chd_gene_set()` with your actual gene data source
2. **Add More Streamers**: Implement streamers for other data sources (gnomAD, expression data, etc.)
3. **Performance Optimization**: Fine-tune chunk sizes based on your datasets
4. **Integration Testing**: Test with real Clinvar data files
5. **Documentation**: Add specific examples for your use cases
