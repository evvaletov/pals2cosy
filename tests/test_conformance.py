"""Conformance with the official PALS examples at the pinned specification commit."""

import builtins
import os
import re
import warnings

import pytest
import yaml

import pals2cosy
from pals2cosy.converter import convert
from pals2cosy.pals_parser import parse_lattice
from tests.cosy_harness import requires_cosy, run_cosy

EXAMPLES = os.path.join(os.path.dirname(__file__), "pals_examples")

with open(os.path.join(EXAMPLES, "manifest.yaml")) as _f:
    MANIFEST = yaml.safe_load(_f)["examples"]

CONVERTING = sorted(k for k, v in MANIFEST.items() if v["expect"] == "converts")
FAILING = sorted(k for k, v in MANIFEST.items() if v["expect"] == "error")


def _vendored_yaml_files():
    found = set()
    for root, _, files in os.walk(EXAMPLES):
        for name in files:
            if name.endswith(".yaml") and name != "manifest.yaml":
                found.add(os.path.relpath(os.path.join(root, name), EXAMPLES))
    return found


def _convert(rel):
    entry = MANIFEST[rel]
    bp, elements = parse_lattice(os.path.join(EXAMPLES, rel), **entry.get("options", {}))
    return bp, elements, convert(bp, elements)


def test_spec_commit_pinned():
    with open(os.path.join(EXAMPLES, "SPEC_COMMIT")) as f:
        commit = f.read().strip()
    assert commit.startswith(pals2cosy.__pals_commit__)
    assert commit == pals2cosy._version.PALS_COMMIT_FULL


def test_manifest_covers_every_example():
    assert _vendored_yaml_files() == set(MANIFEST)
    for rel, entry in MANIFEST.items():
        assert entry["expect"] in ("converts", "error", "fragment"), rel
        if entry["expect"] == "error":
            assert {"error", "match", "reason"} <= set(entry), rel
        if entry["expect"] == "fragment":
            assert rel.endswith(".subpals.yaml"), rel


@pytest.mark.parametrize("rel", CONVERTING)
def test_example_converts(rel):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        bp, elements, fox = _convert(rel)
    assert "PROCEDURE LATTICE ;" in fox
    messages = [str(w.message) for w in caught]
    for pattern in MANIFEST[rel].get("warns", []):
        assert any(re.search(pattern, m) for m in messages), (pattern, messages)


@pytest.mark.parametrize("rel", FAILING)
def test_example_fails_with_documented_reason(rel):
    entry = MANIFEST[rel]
    error = getattr(pals2cosy, entry["error"], None) or getattr(builtins, entry["error"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(error, match=entry["match"]):
            _convert(rel)


@requires_cosy
@pytest.mark.parametrize("rel", CONVERTING)
def test_example_runs_in_cosy(rel, cosy_dir):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bp, elements, fox = _convert(rel)
    result = run_cosy(fox, cosy_dir)
    assert not result.errors, result.errors
    assert result.scalars["SPOS"] == pytest.approx(sum(e["length"] for e in elements))


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_drift_quad_bend_values():
    """Spot-check the derived values the drift_quad_bend example documents."""
    bp, elements, fox = _convert("unit_tests/elements/drift_quad_bend.pals.yaml")
    q1, b1 = elements[1], elements[2]
    assert q1["bn"] == {1: 1.2}
    assert b1["length"] == 1.5 and b1["angle"] == pytest.approx(0.15 * 180 / 3.141592653589793)
    mq = next(l for l in fox.splitlines() if l.strip().startswith("MQ"))
    assert float(mq.split()[2]) == pytest.approx(1.2 * 0.0135)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_iota_ring_closes():
    """The IOTA ring (sector bends, mirrored half via repeat: -1) bends by 2 pi."""
    _, elements, _ = _convert("iota.pals.yaml")
    total = sum(e["angle"] for e in elements if e["type"] == "DPH")
    assert total == pytest.approx(360.0, rel=1e-6)
    names = [e["name"] for e in elements]
    half = names[1:names.index("qe3")]
    assert names[names.index("qe3") + 1:-1] == half[::-1]
