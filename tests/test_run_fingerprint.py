"""Smoke test for running ``ngs-chew fingerprint``"""

import numpy as np

from chew.cli import cli


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
