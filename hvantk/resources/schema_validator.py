"""
Schema validation utilities for hvantk registry datasets.
"""
import json
import jsonschema
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

class SchemaValidator:
    """Validates dataset metadata against appropriate schemas."""

    def __init__(self, schema_dir: Optional[Path] = None):
        self.schema_dir = schema_dir or Path(__file__).parent.parent / "resources" / "schemas"
        self.schemas = {}
        self._load_schemas()

    def _load_schemas(self):
        """Load all available schemas."""
        try:
            for schema_file in self.schema_dir.glob("*.json"):
                schema_name = schema_file.stem
                with open(schema_file, 'r') as f:
                    self.schemas[schema_name] = json.load(f)
                logger.info(f"Loaded schema: {schema_name}")
        except Exception as e:
            logger.error(f"Error loading schemas: {e}")
            raise

    def validate_dataset(self, dataset_metadata: Dict, schema_type: str) -> Tuple[bool, List[str]]:
        """
        Validate dataset metadata against specified schema.

        Args:
            dataset_metadata: Dataset metadata dictionary
            schema_type: Type of schema to validate against

        Returns:
            Tuple of (is_valid, error_messages)
        """
        if schema_type not in self.schemas:
            return False, [f"Schema '{schema_type}' not found. Available: {list(self.schemas.keys())}"]

        try:
            jsonschema.validate(dataset_metadata, self.schemas[schema_type])
            return True, []
        except jsonschema.ValidationError as e:
            return False, [str(e)]
        except Exception as e:
            return False, [f"Validation error: {str(e)}"]

    def get_schema_requirements(self, schema_type: str) -> Dict:
        """Get required fields for a specific schema type."""
        if schema_type not in self.schemas:
            return {}

        schema = self.schemas[schema_type]
        required_fields = schema.get("required", [])

        # Handle allOf schemas (which extend common schema)
        if "allOf" in schema:
            for subschema in schema["allOf"]:
                if "required" in subschema:
                    required_fields.extend(subschema["required"])

        return {
            "required_fields": required_fields,
            "properties": schema.get("properties", {})
        }

    def suggest_schema_type(self, dataset_metadata: Dict) -> str:
        """Suggest appropriate schema type based on dataset metadata."""
        # Simple heuristics to suggest schema type
        if any(key in dataset_metadata for key in ["expression_unit", "platform_type"]):
            return "transcriptomics_schema"
        elif any(key in dataset_metadata for key in ["quantification_method", "protein_inference_method"]):
            return "proteomics_schema"
        elif any(key in dataset_metadata for key in ["genome_build", "variant_types"]):
            return "genomics_schema"
        elif any(key in dataset_metadata for key in ["epigenetic_mark", "assay_type"]):
            return "epigenomics_schema"
        else:
            return "common_schema"
