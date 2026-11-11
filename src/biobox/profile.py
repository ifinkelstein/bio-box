"""Host profile: the resource caps of a shared machine, loaded from YAML (card C-01)."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import yaml

_TYPES = {"int": int, "str": str, "bool": bool}


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
    # Card C-07: true where every step and search worker must run under a memory cap. run-step and
    # the search runner then refuse to start when they can neither open a systemd scope nor find a
    # cgroup memory.max above them (inside the image: launch with scripts/pod_run.sh). Optional.
    require_mem_cap: bool = False
    # Card B-36 (DEC-25): the most CPU-minutes the clustering step of a release may be estimated at in one
    # process; only clustering is exempt from `cpu_min_per_process`, and only under this cap. Optional: a
    # profile without it (or with null) gives no exemption, so clustering over `cpu_min_per_process` is refused.
    cluster_cpu_min_max: int | None = None

    def threads_at(self, hour: int, weekday: int = 0) -> int:
        """Thread cap at a local hour (0-23) and weekday (Monday = 0): the day cap on Monday to
        Friday in [day_start_hour, day_end_hour), the night cap otherwise. The night cap is also
        conditional on no other compute running, which the caller checks (run-step)."""
        day = weekday < 5 and self.day_start_hour <= hour < self.day_end_hour
        return self.max_threads_day if day else self.max_threads_night


def load(path: str | Path) -> Profile:
    """Load and check a profile: every key without a default is required, no unknown keys, and
    each value has its field's type (a bool is not an int, an int is not a bool). A field typed
    `<type> | None` also takes null, which is its default (`build-stage` writes the profile back out
    for its shard processes with every field, nulls included)."""
    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path}: profile must be a mapping")
    fields = {f.name: f for f in dataclasses.fields(Profile)}
    for key in data:
        if key not in fields:
            raise ValueError(f"{path}: unknown key {key!r}")
    for key, f in fields.items():
        if key not in data and f.default is dataclasses.MISSING:
            raise ValueError(f"{path}: missing key {key!r}")
    for key, value in data.items():
        name = fields[key].type if isinstance(fields[key].type, str) else fields[key].type.__name__
        if name.endswith(" | None"):
            if value is None:
                continue
            name = name.removesuffix(" | None")
        want = _TYPES[name]
        if not isinstance(value, want) or (want is not bool and isinstance(value, bool)):
            raise ValueError(f"{path}: key {key!r} must be {want.__name__}, got {type(value).__name__}")
    return Profile(**data)
