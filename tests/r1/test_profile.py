"""C-01: host profile loader (ring R1 by directory)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from biobox import profile

SHIPPED = Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml"


def _raw() -> dict:
    return yaml.safe_load(SHIPPED.read_text())


def _write(tmp_path: Path, data: dict) -> Path:
    p = tmp_path / "p.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


def test_shipped_profile_loads():
    p = profile.load(SHIPPED)
    assert p.max_threads_day == 40
    assert p.max_threads_night == 56
    assert p.mem_scope == "200G"
    assert p.tmp == "/tmp/if2453/bio-box"
    assert p.nice_monolithic == 19
    assert p.nice_fanout == 15
    assert p.cpu_min_per_process == 300
    assert p.small_file_floor_bytes == 64 * 1024 * 1024


def test_unknown_key_refused(tmp_path):
    d = _raw()
    d["surprise"] = 1
    with pytest.raises(ValueError, match="surprise"):
        profile.load(_write(tmp_path, d))


def test_missing_key_refused(tmp_path):
    d = _raw()
    del d["nice_fanout"]
    with pytest.raises(ValueError, match="nice_fanout"):
        profile.load(_write(tmp_path, d))


def test_wrong_type_refused(tmp_path):
    d = _raw()
    d["max_threads_day"] = "forty"
    with pytest.raises(ValueError, match="max_threads_day"):
        profile.load(_write(tmp_path, d))


def test_bool_is_not_an_int(tmp_path):
    d = _raw()
    d["max_threads_night"] = True
    with pytest.raises(ValueError, match="max_threads_night"):
        profile.load(_write(tmp_path, d))


@pytest.mark.parametrize("hour,expected", [(7, 56), (8, 40), (19, 40), (20, 56)])
def test_threads_at(hour, expected):
    assert profile.load(SHIPPED).threads_at(hour) == expected
