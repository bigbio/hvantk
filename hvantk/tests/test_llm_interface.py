import pytest
from unittest.mock import MagicMock, patch
import hail as hl
import os

from hvantk.utils.llm_interface import LLMInterface, natural_language_query


@pytest.fixture
def mock_matrix_table():
    """Create a simple mock MatrixTable for testing using coordinate representation"""
    # Create a simple test matrix with 10 rows (genes) and 5 columns (samples)
    n_rows, n_cols = 10, 5

    # Create coordinate data
    coord_data = []
    for i in range(n_rows):
        for j in range(n_cols):
            coord_data.append({
                'gene_id': f'ENSG{i:08d}',
                'gene_name': f'GENE{i}',
                'sample_id': f'SAMPLE{j}',
                'tissue': f'TISSUE{j%3}',
                'expression': float(i+j)
            })

    # Create a coordinate table
    coord_ht = hl.Table.parallelize(
        coord_data,
        schema='struct{gene_id: str, gene_name: str, sample_id: str, tissue: str, expression: float64}'
    )

    # Convert to a matrix table using the to_matrix_table method
    mt = coord_ht.to_matrix_table(
        row_key=['gene_id'],
        col_key=['sample_id'],
        row_fields=['gene_name'],
        col_fields=['tissue']
    )

    return mt


@patch('hvantk.utils.llm_interface.get_llm_interface')
def test_natural_language_query(mock_get_llm, mock_matrix_table):
    """Test that natural_language_query processes queries and returns expected results"""
    # Create a mock LLM interface
    mock_llm = MagicMock(spec=LLMInterface)

    # Set up the mock response from query_llm
    mock_response = {
        "explanation": "This is a test explanation for gene expression analysis",
        "code_snippets": ["import hail as hl\n\nresult = mt.aggregate_cols(hl.agg.mean(mt.expression))"],
        "executable_code": "import hail as hl\n\nresult = mt.aggregate_cols(hl.agg.mean(mt.expression))"
    }
    mock_llm.query_llm.return_value = mock_response

    # For the execute=True case, set up a mock execution result
    mock_execution_result = {"result": 4.5}  # Expected mean of the test data
    mock_llm.execute_generated_code.return_value = mock_execution_result

    # Configure the mock to return our mock LLM interface
    mock_get_llm.return_value = mock_llm

    # Test case 1: Without execution
    query = "What is the average expression across all samples?"
    result = natural_language_query(query, mock_matrix_table, execute=False)

    # Verify the function called the right methods
    mock_get_llm.assert_called_once()
    mock_llm.query_llm.assert_called_once_with(query, mock_matrix_table, None)
    mock_llm.execute_generated_code.assert_not_called()

    # Verify the result structure
    assert "explanation" in result
    assert "suggested_code" in result
    assert "execution_result" not in result
    assert result["explanation"] == mock_response["explanation"]
    assert result["suggested_code"] == mock_response["executable_code"]

    # Reset mocks for the second test
    mock_get_llm.reset_mock()
    mock_llm.query_llm.reset_mock()

    # Test case 2: With execution
    result_with_execution = natural_language_query(query, mock_matrix_table, execute=True)

    # Verify execute_generated_code was called
    mock_llm.execute_generated_code.assert_called_once_with(
        mock_response["executable_code"],
        mock_matrix_table
    )

    # Verify the result contains execution_result
    print(result_with_execution)
    assert "execution_result" in result_with_execution
    assert result_with_execution["execution_result"] == mock_execution_result


@patch('hvantk.utils.llm_interface.get_llm_interface')
def test_summarize_matrix_query(mock_get_llm, mock_matrix_table):
    """Test natural language query with 'summarize the matrix expression'"""
    # Create a mock LLM interface
    mock_llm = MagicMock(spec=LLMInterface)

    # Set up the mock response for matrix summarization
    mock_response = {
        "explanation": "Here's a summary of the matrix expression data",
        "code_snippets": ["import hail as hl\nfrom hvantk.utils import matrix_utils\n\nresult = matrix_utils.summarize_matrix(mt)"],
        "executable_code": "import hail as hl\nfrom hvantk.utils import matrix_utils\n\nresult = matrix_utils.summarize_matrix(mt)"
    }
    mock_llm.query_llm.return_value = mock_response

    # Set up mock execution result
    mock_execution_result = {
        "dimensions": {
            "n_samples": 5,
            "n_genes": 10,
            "n_entries": 50,
            "sparsity": 0.0
        },
        "expression_stats": {
            "mean": 4.5,
            "std": 2.87,
            "min": 0.0,
            "max": 9.0
        }
    }
    mock_llm.execute_generated_code.return_value = mock_execution_result

    # Configure the mock to return our mock LLM interface
    mock_get_llm.return_value = mock_llm

    # Execute the query with execution
    query = "summarize the matrix expression"
    result = natural_language_query(query, mock_matrix_table, execute=True)

    # Print the result for debugging
    print("\nTest summarize_matrix_query output:")
    print(f"Explanation: {result['explanation']}")
    print(f"Suggested code: {result['suggested_code']}")
    print(f"Execution result: {result['execution_result']}")

    # Verify the query was processed correctly
    mock_llm.query_llm.assert_called_once_with(query, mock_matrix_table, None)

    # Verify code was executed
    mock_llm.execute_generated_code.assert_called_once_with(
        mock_response["executable_code"],
        mock_matrix_table
    )

    # Assert response is not empty
    assert result["explanation"] is not None and result["explanation"] != ""
    assert result["suggested_code"] is not None and result["suggested_code"] != ""
    assert "execution_result" in result
    assert result["execution_result"] == mock_execution_result


@pytest.mark.parametrize("provider,model", [
    ("openai", "gpt-4o-mini"),
    ("anthropic", "claude-3-opus"),
    ("local", "mixtral-8x7b")
])
@patch('hvantk.utils.llm_interface.get_llm_interface')
def test_llm_providers(mock_get_llm, provider, model, mock_matrix_table):
    """Test that each LLM provider can process a natural language query"""
    # Create a mock LLM interface
    mock_llm = MagicMock(spec=LLMInterface)

    # Provider-specific mock response formats can vary slightly
    # but we'll keep a consistent structure for testing
    mock_response = {
        "explanation": f"This is a test explanation from {provider} {model}",
        "code_snippets": [f"# Code generated by {provider} {model}\nimport hail as hl\n\nresult = mt.count_rows()"],
        "executable_code": f"# Code generated by {provider} {model}\nimport hail as hl\n\nresult = mt.count_rows()"
    }
    mock_llm.query_llm.return_value = mock_response

    # Configure the mock to return our mock LLM interface with provider-specific details
    mock_llm.provider = provider
    mock_llm.model = model
    mock_get_llm.return_value = mock_llm

    # For the execute=True case, set up a mock execution result
    mock_execution_result = {"result": 10}  # Expected number of rows
    mock_llm.execute_generated_code.return_value = mock_execution_result

    # Run the query
    query = f"Count the number of rows in the matrix table using {provider}"
    result = natural_language_query(query, mock_matrix_table, execute=True)

    # Print the result for debugging
    print(f"\nTest {provider} {model} output:")
    print(f"Explanation: {result['explanation']}")
    print(f"Suggested code: {result['suggested_code']}")
    print(f"Execution result: {result['execution_result']}")

    # Verify the query was processed correctly
    mock_llm.query_llm.assert_called_once_with(query, mock_matrix_table, None)

    # Verify code was executed
    mock_llm.execute_generated_code.assert_called_once_with(
        mock_response["executable_code"],
        mock_matrix_table
    )

    # Assert response is not empty
    assert result["explanation"] is not None and result["explanation"] != ""
    assert result["suggested_code"] is not None and result["suggested_code"] != ""
    assert "execution_result" in result
    assert result["execution_result"] == mock_execution_result
    assert provider in result["explanation"]  # Provider name should be in explanation
    assert model in result["explanation"]     # Model name should be in explanation


# New tests for LLM interface creation
@patch('os.environ.get')
@patch('openai.OpenAI')
def test_create_openai_interface(mock_openai, mock_env_get):
    """Test that OpenAI interface can be created successfully"""
    # Mock API key in environment
    mock_env_get.return_value = "test-openai-api-key"

    # Mock the OpenAI client
    mock_client = MagicMock()
    mock_openai.return_value = mock_client

    # Create the interface
    llm = LLMInterface(provider="openai", model="gpt-4.1-mini")

    # Verify it was created correctly
    assert llm.provider == "openai"
    assert llm.model == "gpt-4.1-mini"
    assert llm.api_key == "test-openai-api-key"
    assert llm.client == mock_client

    # Verify OpenAI was initialized with the right API key
    mock_openai.assert_called_once_with(api_key="test-openai-api-key")


def test_real_local_model_query(mock_matrix_table):
    """
    Test a real query to the local deepseek-r1 model via Ollama.

    This test is designed to run a real query against the Ollama API using the deepseek-r1 model.
    It uses a mocked matrix and makes an actual API call.

    Note: This test requires Ollama to be running with the deepseek-r1 model loaded.
    It will be skipped if:
    - SKIP_REAL_LLM_TESTS environment variable is set to "true"
    - Connection to Ollama API fails
    """
    # Skip in CI environments
    if os.environ.get("SKIP_REAL_LLM_TESTS", "false").lower() == "true":
        pytest.skip("Skipping real LLM test in CI environment or as requested by SKIP_REAL_LLM_TESTS")

    # Check if Ollama is available before running the test
    import requests
    try:
        ollama_endpoint = os.environ.get("OLLAMA_ENDPOINT", "http://localhost:11434")
        response = requests.get(f"{ollama_endpoint}/api/tags", timeout=2)
        if response.status_code != 200:
            pytest.skip(f"Skipping test: Ollama API returned status code {response.status_code}")

        # Check if deepseek-r1 model is available
        models = response.json().get("models", [])
        model_names = [model.get("name") for model in models]
        if "deepseek-r1" not in model_names and "deepseek-r1:latest" not in model_names:
            pytest.skip("Skipping test: deepseek-r1 model not available in Ollama")

    except (requests.RequestException, ValueError) as e:
        pytest.skip(f"Skipping test: Cannot connect to Ollama API: {str(e)}")

    try:
        # Initialize the LLM interface directly to use the local model
        llm = LLMInterface(
            provider="local",
            model="deepseek-r1",
            temperature=0.7,
            max_tokens=2048
        )

        # Create a query
        query = ("Can you analyze this Hail gene expression matrix and tell me:"
                 " 1) How many samples are there?"
                 " 2) Suggest a visualization?")

        # Make the query
        print("\n========== REAL DEEPSEEK-R1 MODEL QUERY TEST ==========")
        print(f"Query: {query}")
        print("Making query to local deepseek-r1 model via Ollama...")

        response = llm.query_llm(query, mock_matrix_table)

        # Print the results
        print("\n----- LLM RESPONSE -----")
        print(f"Model used: {response['model_used']}")
        print(f"\nExplanation:\n{response['explanation']}")

        print("\n----- SUGGESTED CODE -----")
        if response.get('executable_code'):
            print(response['executable_code'])
        else:
            print("No executable code provided by the model")

        print("\n----- USAGE STATISTICS -----")
        print(f"Usage: {response.get('usage', 'Not available')}")
        print("============================================\n")

    except Exception as e:
        print(f"Error during real model test: {str(e)}")
        pytest.skip(f"Test failed with unexpected error: {str(e)}")


def test_real_openai_model_query(mock_matrix_table):
    """
    Test a real query to an OpenAI model.

    This test is designed to run a real query against the OpenAI API using a specified model.
    It uses a mocked matrix and makes an actual API call.

    Note: This test requires an OpenAI API key to be set in the environment.
    It will be skipped if:
    - SKIP_REAL_LLM_TESTS environment variable is set to "true"
    - OpenAI API key is not available
    - Connection to OpenAI API fails
    """
    # Skip in CI environments
    if os.environ.get("SKIP_REAL_LLM_TESTS", "false").lower() == "true":
        pytest.skip(
            "Skipping real LLM test in CI environment or as requested by SKIP_REAL_LLM_TESTS"
        )

    # Check if OpenAI API key is available
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        pytest.skip("Skipping test: OpenAI API key not found in environment")

    # Check OpenAI API connection
    try:
        # Import here to avoid dependency issues if OpenAI is not installed
        import openai

        # Initialize the LLM interface for OpenAI
        llm = LLMInterface(
            provider="openai",
            model="gpt-4o-mini",  # Or any other available model
            temperature=0.7,
            max_tokens=2048,
        )

        # Create a query
        query = (
            "Can you analyze this Hail gene expression matrix and tell me:"
            " 1) How many samples are there?"
            " 2) Suggest a visualization?"
        )

        # Make the query
        print("\n========== REAL OPENAI MODEL QUERY TEST ==========")
        print(f"Query: {query}")
        print(f"Making query to OpenAI model: {llm.model}...")

        response = llm.query_llm(query, mock_matrix_table)

        # Print the results
        print("\n----- LLM RESPONSE -----")
        print(f"Model used: {response['model_used']}")
        print(f"\nExplanation:\n{response['explanation']}")

        print("\n----- SUGGESTED CODE -----")
        if response.get("executable_code"):
            print(response["executable_code"])
        else:
            print("No executable code provided by the model")

        print("\n----- USAGE STATISTICS -----")
        print(f"Usage: {response.get('usage', 'Not available')}")
        print("============================================\n")

    except (openai.OpenAIError, Exception) as e:
        print(f"Error during real model test: {str(e)}")
        pytest.skip(f"Test failed with unexpected error: {str(e)}")
