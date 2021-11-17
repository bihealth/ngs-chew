"""Smoke test for running ``ngs-chew fingerprint``"""

import numpy as np

from chew.__main__ import main

def test_smoke_test_run_fingerprint_fastq(path_tests, tmpdir):
    main(
        [
            "fingerprint",
            "--reference",
            str(path_tests / "data" / "hs37d5.chr1_fragment.fa"),
            "--output-fingerprint",
            str(tmpdir / "out"),
            "fastq",
            "--input",
            str(path_tests / "data" / "igsr.HG00102.TP73.1.fastq.gz"),
            str(path_tests / "data" / "igsr.HG00102.TP73.2.fastq.gz"),
            "--k",
            "10"
            "--bf_size",
            "1G",
            "depth",
            "10"
        ]
    )

    # Check that output path exists and values equal the reference on in tests/data
    output = tmpdir / "out.npz"
    assert output.exists()

    result = np.load(str(output))
    expected = np.load(str(path_tests / "data/igsr.HG00102.TP73.npz"))
    assert np.array_equal(result, expected)
