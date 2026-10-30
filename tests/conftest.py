"""Records every test execution as `test_run` rows when BIOBOX_TEST_RUN_DIR is set.

    BIOBOX_TEST_RUN_DIR=$ATD_ROOT/tests uv run pytest
    python -m biobox.gate S0 --dir $ATD_ROOT/tests

See src/biobox/testrun.py for how rows are named and stored.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path

import pytest

from biobox import testrun

_ROWS: list[dict] = []
_REPO = Path(__file__).resolve().parents[1]
_OUTCOME = {"passed": "pass", "failed": "fail", "skipped": "skip"}


def pytest_configure(config):
    config.addinivalue_line("markers", "covers(*ids): plan test IDs this test covers (data_layer_and_tests.md)")


def _ring(item) -> str:
    for part in Path(str(item.fspath)).parts:
        if len(part) == 2 and part[0] == "r" and part[1].isdigit():
            return part.upper()
    return "R1"


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    ids = testrun.ids_from_name(item.originalname or item.name)
    for m in item.iter_markers("covers"):
        ids.extend(m.args)
    report.biobox_ids = sorted(set(ids)) or [item.nodeid]
    report.biobox_ring = _ring(item)


def pytest_runtest_logreport(report):
    if not os.environ.get("BIOBOX_TEST_RUN_DIR"):
        return
    if report.when == "call" or (report.when == "setup" and report.outcome != "passed"):
        outcome = _OUTCOME[report.outcome]
        if report.when == "setup" and report.outcome == "failed":
            outcome = "error"
        started = datetime.fromtimestamp(getattr(report, "start", datetime.now().timestamp()), tz=timezone.utc)
        metrics = {k: v for k, v in report.user_properties}
        commit = _commit()
        for tid in getattr(report, "biobox_ids", [report.nodeid]):
            _ROWS.append(testrun.make_row(tid, getattr(report, "biobox_ring", "R1"), outcome, report.duration,
                                          started, report.nodeid, metrics, commit=commit))


_COMMIT: list[str] = []


def _commit() -> str:
    if not _COMMIT:
        _COMMIT.append(testrun.code_commit(_REPO))
    return _COMMIT[0]


def pytest_sessionfinish(session, exitstatus):
    d = os.environ.get("BIOBOX_TEST_RUN_DIR")
    if d and _ROWS:
        path = testrun.write_session(_ROWS, Path(d))
        session.config.get_terminal_writer().line(f"test_run: {len(_ROWS)} rows -> {path}")
