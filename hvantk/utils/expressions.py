"""
Utils general-propose Hail expression
"""

from typing import Union

import hail as hl


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression selecting bi-allelic sites only,
    accounting for whether the input MT/HT was split.
    :param t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    """
    return ~t.was_split if "was_split" in t.row else (hl.len(t.alleles) == 2)


def split_field_expr(
    t: Union[hl.Table, hl.MatrixTable],
    field_name: str = "gene",
    split_char: str = "\\|",
) -> hl.expr.StringExpression:
    """
    Splits the value of a specified field in a Hail Table or MatrixTable based on a delimiter and
    returns the first element of the split result.

    :param t: The Hail Table or MatrixTable containing the field to be split.
    :param field_name: The name of the field to split. Defaults to "gene".
    :param split_char: The delimiter to use for splitting the field's value. Defaults to "\\|".
    :return: The first element of the split string as a StringExpression.
    :rtype: hl.expr.StringExpression
    """
    return t[field_name].split(split_char)[0] if field_name in t.row else t[field_name]
