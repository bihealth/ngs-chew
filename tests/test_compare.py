from chew import compare


def test_load_fingerprint(path_tests):
    header, fingerprint = compare.load_fingerprint(
        False, str(path_tests / "data/igsr.HG00102.TP73.npz")
    )
    assert header == "HG00102-N1-DNA1-WES1"
    assert fingerprint.shape == (3, 23770)
