from pathlib import Path

import pytest

from active_learning.tools import (
    apply_mutation_string,
    parse_fasta_first_sequence,
    sequence_id_from_mutations,
    sequence_sha256,
)


def test_parse_fasta_selects_record_id(tmp_path: Path):
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">A\nACD\n>B\nGGG\n", encoding="utf-8")

    assert parse_fasta_first_sequence(fasta, record_id="A") == "ACD"
    assert parse_fasta_first_sequence(fasta, record_id="B") == "GGG"


def test_apply_unified_mutation_to_single_sequence():
    sequence = "A" * 25 + "G" + "C"

    mutated = apply_mutation_string(sequence, "HG26H", expected_chain="H")

    assert mutated[25] == "H"
    assert len(mutated) == len(sequence)


def test_apply_semicolon_mutations_and_empty_string():
    assert apply_mutation_string("ACD", "") == "ACD"
    assert apply_mutation_string("ACD", "HA1G;HD3E", expected_chain="H") == "GCE"


def test_apply_mutation_rejects_orig_mismatch():
    with pytest.raises(ValueError, match="HC2A"):
        apply_mutation_string("AAD", "HC2A")


def test_apply_mutation_rejects_unexpected_chain_when_not_ignored():
    with pytest.raises(ValueError, match="LC2A"):
        apply_mutation_string(
            "ACD",
            "LC2A",
            expected_chain="H",
            ignore_chain_for_single_sequence=False,
        )


def test_sequence_identity_helpers_are_stable():
    assert sequence_id_from_mutations("sdab", "HG26H") == "sdab:HG26H"
    assert sequence_id_from_mutations("sdab", "") == "sdab:WT"
    assert len(sequence_sha256("ACD")) == 64
    assert sequence_sha256("ACD") == sequence_sha256("ACD")
