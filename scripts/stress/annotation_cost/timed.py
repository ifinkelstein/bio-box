"""Run one command and record its wall time, CPU time and peak RSS (card 0-09a).

    python timed.py <record.json> -- <command...>

Runs inside the image (the image has no /usr/bin/time). CPU time is the rusage that `os.wait4` returns for
the child, which on Linux includes every descendant the child waited for (MMseqs2's workflow scripts wait
for each module they start). Peak RSS is the rusage maxrss (the largest single process in that tree).
Writes the record atomically and exits with the command's exit code.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import subprocess
import sys
import time


def main(argv: list[str]) -> int:
    if len(argv) < 3 or argv[1] != "--":
        print("usage: timed.py <record.json> -- <command...>", file=sys.stderr)
        return 2
    rec_path, cmd = argv[0], argv[2:]
    started = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    t0 = time.monotonic()
    p = subprocess.Popen(cmd)
    _, status, ru = os.wait4(p.pid, 0)
    wall = time.monotonic() - t0
    code = os.waitstatus_to_exitcode(status)
    rec = {
        "argv": cmd,
        "host": os.uname().nodename,
        "started_at": started,
        "wall_s": round(wall, 3),
        "cpu_min": round((ru.ru_utime + ru.ru_stime) / 60.0, 4),
        "user_min": round(ru.ru_utime / 60.0, 4),
        "sys_min": round(ru.ru_stime / 60.0, 4),
        "peak_rss_mib": round(ru.ru_maxrss / 1024.0, 1),
        "exit_code": code,
    }
    tmp = f"{rec_path}.tmp-{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=1, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, rec_path)
    print(f"timed: {json.dumps(rec, sort_keys=True)}", flush=True)
    return code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
