"""Smoke test for running ``ngs-chew fingerprint``

The unmarked tests use small files such that the test can run quickly.  The tests marked as slow
use large files that taker longer to run and that need to be enabled explicitely.
"""

from pathlib import Path

from _pytest._py.path import LocalPath
from click.testing import CliRunner
import numpy as np
import pytest

from chew.cli import cli
from chew.compare import relatedness


def test_smoke_test_run_fingerprint(cli_runner, path_tests, tmpdir):
    result = cli_runner.invoke(
        cli,
        [
            "fingerprint",
            "--reference",
            str(path_tests / "data" / "hs37d5.chr1_fragment.fa"),
            "--input-bam",
            str(path_tests / "data" / "igsr.HG00102.TP73.bam"),
            "--output-fingerprint",
            str(tmpdir / "out"),
        ],
    )

    # Check that output path exists and values equal the reference on in tests/data
    output = tmpdir / "out.npz"
    assert output.exists()

    result = np.load(str(output))
    expected = np.load(str(path_tests / "data/igsr.HG00102.TP73.npz"))
    assert np.array_equal(result, expected)


@pytest.mark.slow
@pytest.mark.parametrize(
    "path_ref,path_bam",
    [
        (
            "tests/data/fingerprint/hs37d5.fa.gz",
            "tests/data/fingerprint/HG00138_WES_GRCh37.excerpt.bam",
        ),
        (
            "tests/data/fingerprint/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz",
            "tests/data/fingerprint/HG00138_WES_GRCh38.excerpt.bam",
        ),
        (
            "tests/data/fingerprint/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz",
            "tests/data/fingerprint/HG00138_WGS_GRCh38.excerpt.bam",
        ),
    ],
)
def test_fingerprint_bam_grch37(
    cli_runner: CliRunner, tmpdir: LocalPath, path_tests: Path, path_ref: str, path_bam: str
):
    cli_runner.invoke(
        cli,
        [
            "fingerprint",
            "--reference",
            path_ref,
            "--input-bam",
            path_bam,
            "--output-fingerprint",
            str(tmpdir / "out"),
        ],
    )

    # Check that output path exists and is similar to all finger prints.
    output = tmpdir / "out.npz"
    assert output.exists()
    result_np = np.load(str(output))["autosomal_fingerprint"]

    for fname in (
        "HG00138_WES_GRCh37.fingerprint.npz",
        "HG00138_WES_GRCh38.fingerprint.npz",
        "HG00138_WGS_GRCh38.fingerprint.npz",
    ):
        expected_np = np.load(str(path_tests / "data/fingerprint" / fname))["autosomal_fingerprint"]
        n_ibs0, rel, _ = relatedness(result_np, expected_np)
        assert n_ibs0 <= 10
        assert rel > 0.98
