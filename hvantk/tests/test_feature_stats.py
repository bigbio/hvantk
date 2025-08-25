import hail as hl
import pytest

def safe_float(val):
    return hl.missing(hl.tfloat64) if val is None else float(val)

def test_feature_coverage_stats():
    hl.init(log='/tmp/hail_test.log', quiet=True)
    rows = [
        {'combined_deleteriousness': 1.2, 'median_expression': 5.0, 'constraint_score': 0.8, 'AF': 0.01, 'feature_completeness': 0.9},
        {'combined_deleteriousness': None, 'median_expression': 2.0, 'constraint_score': 0.7, 'AF': 0.02, 'feature_completeness': 0.8},
        {'combined_deleteriousness': None, 'median_expression': 3.0, 'constraint_score': None, 'AF': 0.03, 'feature_completeness': 0.7},
        {'combined_deleteriousness': 0.5, 'median_expression': None, 'constraint_score': None, 'AF': None, 'feature_completeness': 0.6},
        {'combined_deleteriousness': None, 'median_expression': None, 'constraint_score': None, 'AF': None, 'feature_completeness': 0.5},
    ]
    ht = hl.Table.parallelize([
        hl.struct(
            combined_deleteriousness=safe_float(row.get('combined_deleteriousness')),
            median_expression=safe_float(row.get('median_expression')),
            constraint_score=safe_float(row.get('constraint_score')),
            AF=safe_float(row.get('AF')),
            feature_completeness=safe_float(row.get('feature_completeness')),
        ) for row in rows
    ])
    row_fields = set(ht.row.dtype.fields)
    feature_stats = ht.aggregate(
        hl.struct(
            has_prediction_scores=(
                hl.agg.fraction(hl.is_defined(ht.combined_deleteriousness))
                if "combined_deleteriousness" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            has_expression=(
                hl.agg.fraction(hl.is_defined(ht.median_expression))
                if "median_expression" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            has_constraint=(
                hl.agg.fraction(hl.is_defined(ht.constraint_score))
                if "constraint_score" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            has_frequency=(
                hl.agg.fraction(hl.is_defined(ht.AF))
                if "AF" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            avg_feature_completeness=(
                hl.agg.mean(ht.feature_completeness)
                if "feature_completeness" in row_fields
                else hl.agg.mean(hl.missing(hl.tfloat64))
            ),
        )
    )
    print("Feature coverage statistics:")
    print(f"  Prediction scores: {feature_stats.has_prediction_scores:.2%}")
    print(f"  Gene expression: {feature_stats.has_expression:.2%}")
    print(f"  Constraint metrics: {feature_stats.has_constraint:.2%}")
    print(f"  Population frequency: {feature_stats.has_frequency:.2%}")
    print(f"  Average feature completeness: {feature_stats.avg_feature_completeness:.2f}")
    # Assert that the output is as expected (basic checks)
    assert 0 <= feature_stats.has_prediction_scores <= 1
    assert 0 <= feature_stats.has_expression <= 1
    assert 0 <= feature_stats.has_constraint <= 1
    assert 0 <= feature_stats.has_frequency <= 1
    assert 0.0 < feature_stats.avg_feature_completeness <= 1.0

