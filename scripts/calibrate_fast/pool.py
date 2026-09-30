# SPDX-License-Identifier: Apache-2.0
"""The worker pool that runs trials (MEDS_FAST_CALIBRATION_PLAN.md §7 P1 "workers").

A trial is one single-threaded meds_main process, and the fit runs hundreds at a time.

- LocalPool: N concurrent subprocesses on this machine (a workstation, or one node).
- QueuePool: a directory queue shared by workers on many nodes. The driver writes one task file
  per trial; each worker (`calibrate_fast.py worker`, one per node, started inside the same Slurm
  allocation) claims tasks by an atomic rename and runs them on its own cores, writing a done
  file. One allocation holds the workers for the whole fit, so no trial waits in the Slurm queue,
  and no per-trial job step loads the scheduler.

Both run a batch of tasks and return {task id: (status, seconds)}; status is "ok", "timeout" or
"exit <code>". A task is (id, argv, cwd, log path, timeout seconds).
"""
from __future__ import annotations

import json
import os
import socket
import subprocess
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Task:
    id: str
    argv: list
    cwd: str
    log: str
    timeout: float


def run_task(t: Task) -> tuple[str, float]:
    t0 = time.time()
    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        with open(t.log, "w") as fh:
            res = subprocess.run(t.argv, cwd=t.cwd, stdout=fh, stderr=subprocess.STDOUT,
                                 timeout=t.timeout, env=env)
        status = "ok" if res.returncode == 0 else f"exit {res.returncode}"
    except subprocess.TimeoutExpired:
        status = "timeout"
    return status, time.time() - t0


class LocalPool:
    """One executor for the pool's life, so concurrent callers (the fit's parallel starts) share
    its workers instead of multiplying them."""

    def __init__(self, workers: int):
        self.workers = max(1, int(workers))
        self.ex = ThreadPoolExecutor(self.workers)

    def run(self, tasks: list[Task]) -> dict:
        futures = [self.ex.submit(run_task, t) for t in tasks]
        return {t.id: f.result() for t, f in zip(tasks, futures)}

    def close(self):
        self.ex.shutdown(wait=True)


class QueuePool:
    """The driver's side of the directory queue at `root` (queue/, claimed/, done/)."""

    def __init__(self, root, poll: float = 0.5):
        self.root = Path(root)
        for sub in ("queue", "claimed", "done"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        self.poll = poll

    def run(self, tasks: list[Task]) -> dict:
        batch = uuid.uuid4().hex[:8]
        ids = {}
        for t in tasks:
            name = f"{batch}-{len(ids):06d}"
            ids[name] = t.id
            tmp = self.root / "queue" / f".{name}.tmp"
            tmp.write_text(json.dumps(t.__dict__))
            tmp.rename(self.root / "queue" / f"{name}.json")
        out = {}
        while len(out) < len(ids):
            for name in list(ids):
                if ids[name] in out:
                    continue
                done = self.root / "done" / f"{name}.json"
                if done.exists():
                    try:
                        rec = json.loads(done.read_text())
                    except json.JSONDecodeError:
                        continue                     # still being written
                    out[ids[name]] = (rec["status"], rec["seconds"])
                    done.unlink()
            if len(out) < len(ids):
                time.sleep(self.poll)
        return out

    def close(self):
        (self.root / "STOP").write_text("stop\n")


def worker(root, slots: int, poll: float = 0.5) -> None:
    """A node's worker: claim tasks from root/queue until root/STOP exists."""
    root = Path(root)
    me = f"{socket.gethostname()}-{os.getpid()}"
    for sub in ("queue", "claimed", "done"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    running = {}
    with ThreadPoolExecutor(slots) as ex:
        while True:
            for name, fut in list(running.items()):
                if fut.done():
                    status, secs = fut.result()
                    tmp = root / "done" / f".{name}.tmp"
                    tmp.write_text(json.dumps({"status": status, "seconds": secs, "worker": me}))
                    tmp.rename(root / "done" / f"{name}.json")
                    (root / "claimed" / f"{me}.{name}.json").unlink(missing_ok=True)
                    del running[name]
            free = slots - len(running)
            if free > 0:
                for q in sorted((root / "queue").glob("*.json"))[:free]:
                    claim = root / "claimed" / f"{me}.{q.name}"
                    try:
                        q.rename(claim)              # atomic: one worker wins
                    except (FileNotFoundError, OSError):
                        continue
                    t = Task(**json.loads(claim.read_text()))
                    running[q.stem] = ex.submit(run_task, t)
            if (root / "STOP").exists() and not running and not any((root / "queue").glob("*.json")):
                return
            time.sleep(poll)


def make_pool(kind: str, workers: int, root=None):
    if kind == "local":
        return LocalPool(workers)
    if kind == "queue":
        return QueuePool(root)
    raise ValueError(f"unknown pool '{kind}'")
