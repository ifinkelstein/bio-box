"""Host profile: the resource caps of a shared machine, loaded from YAML (card C-01)."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import yaml


@dataclasses.dataclass(frozen=True)
class Profile:
    max_threads_day: int
    max_threads_night: int
    mem_scope: str
    tmp: str
    durable: str
    bulk: str
    nice_monolithic: int
    nice_fanout: int
    day_start_hour: int
    day_end_hour: int
    cpu_min_per_process: int
    small_file_floor_bytes: int

    def threads_at(self, hour: int, weekday: int = 0) -> int:
        """Thread cap at a local hour (0-23) and weekday (Monday = 0): the day cap on Monday to
        Friday in [day_start_hour, day_end_hour), the night cap otherwise. The night cap is also
        conditional on no other compute running, which the caller checks (run-step)."""
        day = weekday < 5 and self.day_start_hour <= hour < self.day_end_hour
        return self.max_threads_day if day else self.max_threads_night


def load(path: str | Path) -> Profile:
    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path}: profile must be a mapping")
    fields = {f.name: f.type for f in dataclasses.fields(Profile)}
    for key in data:
        if key not in fields:
            raise ValueError(f"{path}: unknown key {key!r}")
    for key in fields:
        if key not in data:
            raise ValueError(f"{path}: missing key {key!r}")
    for key, annotation in fields.items():
        want = int if annotation in (int, "int") else str
        value = data[key]
        if not isinstance(value, want) or isinstance(value, bool):
            raise ValueError(f"{path}: key {key!r} must be {want.__name__}, got {type(value).__name__}")
    return Profile(**data)
