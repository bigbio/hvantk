"""
LLM Interface for Hail MatrixTable gene expression queries.

This module provides functionality for:
- Setting up connections to various LLM providers (OpenAI, Anthropic, etc.)
- Processing MatrixTable data to be suitable for LLM consumption
- Converting natural language queries into Hail operations
- Handling responses from LLMs and converting them back into executable code

The module is designed to work with the MatrixTable structure defined in matrix_utils.py
and provides a bridge between natural language and bioinformatic analysis.
"""

import json
import logging
import os
from typing import Dict, Optional

import hail as hl
import numpy as np
import pandas as pd

# Import our own utilities
from ..utils import matrix_utils

# Setup basic logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# LLM provider options
LLM_PROVIDERS = ["openai", "anthropic", "local", "google"]


class LLMConfigError(Exception):
    """Raised when there is an issue with the LLM configuration."""

    pass


class LLMInterface:
    """Interface for interacting with Language Models for MatrixTable analysis."""

    def __init__(
        self,
        provider: str = "google",
        model: str = None,
        api_key: str = None,
        max_tokens: int = 4096,
        temperature: float = 0.7,
        matrix_sample_size: int = 1000,
    ):
        """
        Initialize the LLM interface.

        Args:
            provider: LLM provider name (google, openai, anthropic, local)
            model: Specific model to use (defaults to provider's recommended)
            api_key: API key for the provider (if None, looks for env variables)
            max_tokens: Maximum tokens to use in prompts
            temperature: Temperature setting for response generation
            matrix_sample_size: Default number of rows/cols to sample from matrices
        """
        self.provider = provider.lower()
        if self.provider not in LLM_PROVIDERS:
            raise LLMConfigError(
                f"Unsupported LLM provider: {provider}. Choose from {LLM_PROVIDERS}"
            )

        self.model = model or self._get_default_model()
        self.api_key = api_key or self._get_api_key()
        self.max_tokens = max_tokens
        self.temperature = temperature
        self.matrix_sample_size = matrix_sample_size
        self.client = self._initialize_client()

    def _get_default_model(self) -> str:
        """Get the default model for the selected provider."""
        defaults = {
            "openai": "gpt-4o-mini",
            "anthropic": "claude-3-opus",
            "local": "deepseek-r1",
            "google": "gemini-2.0-flash",
        }
        return defaults.get(self.provider, "gemini-2.0-flash")

    def _get_api_key(self) -> Optional[str]:
        """Get API key from environment variables."""
        if self.provider == "openai":
            return os.environ.get("OPENAI_API_KEY")
        elif self.provider == "anthropic":
            return os.environ.get("ANTHROPIC_API_KEY")
        elif self.provider == "google":
            return os.environ.get("GEMINI_API_KEY")
        elif self.provider == "local":
            # Local models typically don't need an API key
            return None

        return None

    def _initialize_client(self):
        """Initialize the appropriate client for the selected LLM provider."""
        if self.provider == "openai":
            try:
                import openai

                client = openai.OpenAI(api_key=self.api_key)
                return client
            except ImportError as e:
                raise LLMConfigError(
                    "OpenAI package not installed. Run 'pip install openai'"
                ) from e
        elif self.provider == "anthropic":
            try:
                import anthropic

                client = anthropic.Anthropic(api_key=self.api_key)
                return client
            except ImportError as e:
                raise LLMConfigError(
                    "Anthropic package not installed. Run 'pip install anthropic'"
                ) from e
        elif self.provider == "google":
            try:
                from google import genai

                genai.configure(api_key=self.api_key)
                return genai
            except ImportError as e:
                raise LLMConfigError(
                    "Google Generative AI package not installed. Run 'pip install google-generativeai'"
                ) from e
        elif self.provider == "local":
            try:
                import requests

                # For Ollama models, we don't need to load the model directly
                # We'll just verify we can connect to the Ollama API
                # Default Ollama endpoint
                ollama_endpoint = os.environ.get(
                    "OLLAMA_ENDPOINT", "http://localhost:11434"
                )
                # Test connection to Ollama
                try:
                    response = requests.get(f"{ollama_endpoint}/api/tags")
                    if response.status_code != 200:
                        raise LLMConfigError(
                            f"Could not connect to Ollama API at {ollama_endpoint}"
                        )
                    return {"endpoint": ollama_endpoint}
                except requests.exceptions.RequestException as e:
                    raise LLMConfigError(
                        f"Failed to connect to Ollama API: {str(e)}"
                    ) from e
            except ImportError as e:
                raise LLMConfigError(
                    "Requests package not installed. Run 'pip install requests'"
                ) from e

        raise LLMConfigError(f"Unsupported LLM provider: {self.provider}")

    def prepare_matrix_data(
        self, mt: hl.MatrixTable, sample_rows: int = None, sample_cols: int = None
    ) -> Dict:
        """
        Prepare MatrixTable data for LLM consumption by sampling and summarizing.

        Args:
            mt: Hail MatrixTable to prepare
            sample_rows: Number of rows to sample (None = use default)
            sample_cols: Number of columns to sample (None = use default)

        Returns:
            Dict containing summarized and sampled data suitable for LLM
        """
        sample_rows = sample_rows or self.matrix_sample_size
        sample_cols = sample_cols or self.matrix_sample_size

        # Get basic summary statistics
        summary = matrix_utils.summarize_matrix(mt)

        # Sample a subset of data for the LLM
        n_rows = mt.count_rows()
        n_cols = mt.count_cols()

        row_fraction = min(1.0, sample_rows / n_rows) if n_rows > 0 else 0
        col_fraction = min(1.0, sample_cols / n_cols) if n_cols > 0 else 0

        # Sample the matrix
        sampled_mt = mt
        if row_fraction < 1.0:
            sampled_mt = sampled_mt.sample_rows(row_fraction)
        if col_fraction < 1.0:
            sampled_mt = sampled_mt.sample_cols(col_fraction)

        # Convert a small sample to pandas for easy display
        # This is a simplified approach - real implementation might be more sophisticated
        sampled_data = sampled_mt.entries().head(100)
        if sampled_data:
            sample_df = sampled_data.to_pandas()
            # Convert to dictionary representation for JSON serialization
            sample_dict = sample_df.to_dict(orient="records")
        else:
            sample_dict = []

        # Extract column and row metadata for context
        col_meta = sampled_mt.col.collect()[:10]  # First 10 samples
        row_meta = sampled_mt.row.collect()[:10]  # First 10 genes

        prepared_data = {
            "summary": summary,
            "col_metadata_sample": [dict(x) for x in col_meta],
            "row_metadata_sample": [dict(x) for x in row_meta],
            "data_sample": sample_dict,
            "matrix_schema": str(mt.describe()),
        }

        return prepared_data

    def _query_openai(self, system_prompt: str, user_prompt: str) -> Dict:
        """Query OpenAI API."""
        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=self.temperature,
            max_tokens=self.max_tokens,
        )
        return {
            "response_text": response.choices[0].message.content,
            "model_used": self.model,
            "finish_reason": response.choices[0].finish_reason,
            "usage": (
                response.usage.dict()
                if hasattr(response.usage, "dict")
                else vars(response.usage)
            ),
        }

    def _query_anthropic(self, system_prompt: str, user_prompt: str) -> Dict:
        """Query Anthropic API."""
        response = self.client.messages.create(
            model=self.model,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}],
            temperature=self.temperature,
            max_tokens=self.max_tokens,
        )
        return {
            "response_text": response.content[0].text,
            "model_used": self.model,
            "stop_reason": response.stop_reason,
            "usage": response.usage,
        }

    def _query_google(self, system_prompt: str, user_prompt: str) -> Dict:
        """Query Google Gemini API."""
        generation_config = {
            "temperature": self.temperature,
            "max_output_tokens": self.max_tokens,
            "systemInstruction": {"parts": [system_prompt]},
        }

        response = self.client.GenerativeModel(
            model_name=self.model, generation_config=generation_config
        ).generate_content(contents=[{"role": "user", "parts": [user_prompt]}])

        return {
            "response_text": response.text,
            "model_used": self.model,
            "finish_reason": "stop",  # Gemini doesn't provide this explicitly
            "usage": {
                "prompt_tokens": -1,
                "completion_tokens": -1,
                "total_tokens": -1,
            },  # Not provided by Gemini API
        }

    def _query_local(self, system_prompt: str, user_prompt: str) -> Dict:
        """Query local LLM (Ollama) API."""
        try:
            import requests

            # Get Ollama client configuration
            ollama_endpoint = self.client["endpoint"]

            # Prepare the request payload for Ollama API
            payload = {
                "model": self.model,
                "prompt": f"{system_prompt}\n\n{user_prompt}",
                "stream": False,
                "options": {
                    "temperature": self.temperature,
                    "num_predict": self.max_tokens,
                },
            }

            # Make the API call to Ollama
            response = requests.post(f"{ollama_endpoint}/api/generate", json=payload)

            if response.status_code != 200:
                raise Exception(f"Ollama API error: {response.text}")

            response_data = response.json()
            response_text = response_data.get("response", "")

            return {
                "response_text": response_text,
                "model_used": self.model,
                "finish_reason": "stop",
                "usage": {
                    "prompt_tokens": response_data.get("prompt_eval_count", -1),
                    "completion_tokens": response_data.get("eval_count", -1),
                    "total_tokens": response_data.get("prompt_eval_count", 0)
                    + response_data.get("eval_count", 0),
                },
            }
        except Exception as e:
            logger.error(f"Error querying local model: {str(e)}")
            return {
                "response_text": "",
                "model_used": self.model,
                "finish_reason": "error",
                "usage": {
                    "prompt_tokens": -1,
                    "completion_tokens": -1,
                    "total_tokens": -1,
                },
                "error": str(e),
            }

    def query_llm(
        self,
        natural_language_query: str,
        mt: hl.MatrixTable,
        context: Optional[Dict] = None,
    ) -> Dict:
        """
        Send a natural language query about a MatrixTable to the LLM.
        Args:
            natural_language_query: The query in natural language
            mt: The MatrixTable to analyze
            context: Additional context to provide to the LLM
        Returns:
            Dict containing LLM response and suggested code
        """
        # Prepare matrix data
        matrix_data = self.prepare_matrix_data(mt)

        # Construct prompt
        system_prompt = self._get_system_prompt()
        user_prompt = self._construct_user_prompt(
            natural_language_query, matrix_data, context
        )

        # Call appropriate LLM based on provider
        provider_methods = {
            "openai": self._query_openai,
            "anthropic": self._query_anthropic,
            "google": self._query_google,
            "local": self._query_local,
        }

        if self.provider not in provider_methods:
            raise LLMConfigError(f"Unsupported LLM provider: {self.provider}")

        result = provider_methods[self.provider](system_prompt, user_prompt)

        # Parse the response to extract code snippets and explanation
        parsed_response = self._parse_llm_response(result["response_text"])
        result.update(parsed_response)

        return result

    def _get_system_prompt(self) -> str:
        """Generate the system prompt for the LLM."""
        return """You are a bioinformatics expert specializing in gene expression analysis using Hail.
                  Your task is to help users analyze gene expression data stored in Hail MatrixTable format.
                  When given a query, respond with:
                  1. A clear explanation of how to solve the problem
                  2. Executable Python code using Hail that implements the solution
                  3. Any relevant insights about the biological implications

                  Use the provided MatrixTable schema and sample data to understand the structure.
                  Always focus on practical, executable solutions with proper Hail syntax.
                  For complex analyses, break down the approach into clear steps.
        """

    def _construct_user_prompt(
        self, query: str, matrix_data: Dict, context: Optional[Dict] = None
    ) -> str:
        """
        Construct the user prompt from query and context.

        Args:
            query: Natural language query
            matrix_data: Prepared matrix data
            context: Additional context

        Returns:
            Formatted prompt text
        """
        context_str = json.dumps(context, default=str) if context else "{}"

        # Create a condensed version of matrix_data to avoid token limit issues
        condensed_data = {
            "summary": matrix_data["summary"],
            "schema": matrix_data["matrix_schema"],
            "sample_rows": matrix_data["row_metadata_sample"][
                :3
            ],  # Limit samples to save tokens
            "sample_cols": matrix_data["col_metadata_sample"][:3],
        }

        matrix_str = json.dumps(condensed_data, default=str, indent=2)

        prompt = f"""
                  I want to analyze a gene expression dataset stored as a Hail MatrixTable.

                  QUERY:
                  {query}

                  MATRIX INFORMATION:
                  {matrix_str}

                  ADDITIONAL CONTEXT:
                  {context_str}

                  Please provide:
                  1. An explanation of how to approach this analysis
                  2. Executable Python code using Hail
                  3. Any insights about the biological implications
        """
        return prompt

    def _parse_llm_response(self, response_text: str) -> Dict:
        """
        Parse the LLM response to extract code snippets and explanation.

        Args:
            response_text: Raw text response from LLM

        Returns:
            Dict with parsed components
        """
        import re

        # Extract code blocks
        code_blocks = re.findall(
            r"```(?:python)?\s*(.*?)\s*```", response_text, re.DOTALL
        )

        # Extract explanation (text outside code blocks)
        explanation = response_text
        for block in code_blocks:
            explanation = explanation.replace(f"```python\n{block}\n```", "")
            explanation = explanation.replace(f"```\n{block}\n```", "")

        # Clean up explanation
        explanation = re.sub(r"\n{3,}", "\n\n", explanation).strip()

        return {
            "explanation": explanation,
            "code_snippets": code_blocks,
            "executable_code": "\n\n".join(code_blocks) if code_blocks else None,
        }

    def execute_generated_code(self, code, mt=None, globals=None):
        """
        Safely execute the code generated by the LLM using RestrictedPython sandbox.

        Args:
            code: Python code to execute
            mt: MatrixTable to operate on
            globals: Additional global variables to provide

        Returns:
            Result of the code execution
        """
        try:
            from RestrictedPython import compile_restricted, safe_globals
            from RestrictedPython.Guards import (
                guarded_getattr,
                guarded_getitem,
                guarded_setitem,
                guarded_iter,
            )
            from RestrictedPython.PrintCollector import PrintCollector
            import time
        except ImportError:
            logger.error(
                "RestrictedPython not installed. Run 'pip install RestrictedPython'"
            )
            return {"error": "RestrictedPython not installed", "code": code}

        # Setup execution environment
        if globals is None:
            globals = {}

        # Create restricted globals
        restricted_globals = safe_globals.copy()

        # Add necessary builtins
        restricted_globals["_print_"] = PrintCollector
        restricted_globals["_getattr_"] = guarded_getattr
        restricted_globals["_getitem_"] = guarded_getitem
        restricted_globals["_setitem_"] = guarded_setitem
        restricted_globals["_iter_"] = guarded_iter

        # Add safe subset of Python builtins
        safe_builtins = {
            "abs": abs,
            "all": all,
            "any": any,
            "bool": bool,
            "dict": dict,
            "enumerate": enumerate,
            "filter": filter,
            "float": float,
            "frozenset": frozenset,
            "int": int,
            "isinstance": isinstance,
            "len": len,
            "list": list,
            "map": map,
            "max": max,
            "min": min,
            "range": range,
            "round": round,
            "set": set,
            "sorted": sorted,
            "str": str,
            "sum": sum,
            "tuple": tuple,
            "zip": zip,
        }
        restricted_globals["__builtins__"].update(safe_builtins)

        # Add specific whitelisted modules and objects
        allowed_modules = {
            "hl": hl,
            "np": np,
            "pd": pd,
            "mt": mt,
            "matrix_utils": matrix_utils,
        }
        restricted_globals.update(allowed_modules)
        restricted_globals.update(globals)

        # Set execution timeout (30 seconds)
        timeout = 30.0
        start_time = time.time()

        # Add a hook to check for timeouts during execution
        def check_timeout():
            if time.time() - start_time > timeout:
                raise TimeoutError(f"Execution timed out (> {timeout} seconds)")
            return True

        restricted_globals["_check_timeout"] = check_timeout

        try:
            # Add timeout checks to the code
            modified_code = "if _check_timeout():\n"
            for line in code.split("\n"):
                modified_code += f"    {line}\n    if _check_timeout(): pass\n"

            # Compile the code in restricted mode
            byte_code = compile_restricted(
                modified_code, filename="<generated_code>", mode="exec"
            )

            # Create a namespace for execution
            local_vars = {}

            # Execute the compiled code in the restricted environment
            exec(byte_code, restricted_globals, local_vars)

            # Look for a result variable - common convention in generated code
            if "result" in local_vars:
                return local_vars["result"]
            elif "_print" in local_vars:  # Get any printed output
                return local_vars["_print"]()
            else:
                # Return all local variables as the result (excluding internal ones)
                return {
                    k: v
                    for k, v in local_vars.items()
                    if not k.startswith("_") and k not in restricted_globals
                }

        except TimeoutError as te:
            logger.error(f"Timeout executing generated code: {str(te)}")
            return {"error": str(te), "code": code}
        except SyntaxError as se:
            logger.error(f"Syntax error in generated code: {str(se)}")
            return {"error": f"Syntax error: {str(se)}", "code": code}
        except Exception as e:
            logger.error(f"Error executing generated code: {str(e)}")
            return {"error": str(e), "code": code}


def get_llm_interface(**kwargs) -> LLMInterface:
    """
    Get an LLM interface with default or specified configuration.

    This is a factory function that creates and returns an LLMInterface instance
    with either default settings or the specified configuration.

    Args:
        **kwargs: Configuration options to pass to LLMInterface constructor

    Returns:
        An initialized LLMInterface instance
    """
    return LLMInterface(**kwargs)


def natural_language_query(
    query: str,
    mt: hl.MatrixTable,
    context: Optional[Dict] = None,
    execute: bool = False,
    llm_config: Optional[Dict] = None,
) -> Dict:
    """
    Process a natural language query about a MatrixTable and return results.

    This function provides a simplified interface for querying an LLM about a Hail MatrixTable.
    It handles the setup of the LLM, sends the query, and optionally executes the generated code.

    Args:
        query: Natural language query about the MatrixTable
        mt: The Hail MatrixTable to analyze
        context: Additional context to provide to the LLM (optional)
        execute: Whether to execute the generated code (default: False)
        llm_config: Configuration for the LLM interface (optional)

    Returns:
        Dict containing explanation, suggested code, and optionally execution results
    """
    # Initialize the LLM interface with provided or default configuration
    llm_config = llm_config or {}
    llm = get_llm_interface(**llm_config)

    # Send the query to the LLM
    response = llm.query_llm(query, mt, context)

    # Prepare the result dictionary
    result = {
        "explanation": response.get("explanation", ""),
        "suggested_code": response.get("executable_code", ""),
    }

    # Execute the generated code if requested
    if execute and response.get("executable_code"):
        result["execution_result"] = llm.execute_generated_code(
            response["executable_code"], mt
        )

    return result
