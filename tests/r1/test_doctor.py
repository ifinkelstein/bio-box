"""Ring R1: `biobox.doctor` produces a complete, JSON-serializable report."""

import json

from biobox import doctor


def test_report_has_every_section():
    r = doctor.report()
    assert set(r) >= {"python", "platform", "cpu", "memory", "ulimits", "filesystems", "tools", "packages"}
    assert r["packages"]["pyarrow"]
    json.dumps(r, default=str)
