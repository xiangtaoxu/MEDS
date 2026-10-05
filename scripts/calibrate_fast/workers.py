# SPDX-License-Identifier: Apache-2.0
"""The workers that run trials and chains: each is one single-threaded MEDS process, and a fit runs
hundreds at a time.

- LocalWorkers: N processes at once on this machine (a workstation, or one node).
- QueueWorkers: a directory queue shared by workers on many nodes. The driver writes one task file
  per run; each node's `calibrate_fast.py worker` (one per node, started inside the same Slurm
  allocation) claims tasks by an atomic rename, runs them on its own cores and writes a done file.
  One allocation holds the workers for the whole fit, so no run waits in the Slurm queue.

Both run a batch of tasks and return {task id: (status, seconds)}, status "ok", "timeout" or
"exit <code>". A task is (id, command, directory, log file, timeout in seconds).
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


class LocalWorkers:
    """N processes at once on this machine, shared by every caller (the chains start together)."""

    def __init__(self, n: int):
        self.ex = ThreadPoolExecutor(max(1, int(n)))

    def run(self, tasks: list[Task]) -> dict:
        futures = [self.ex.submit(run_task, t) for t in tasks]
        return {t.id: f.result() for t, f in zip(tasks, futures)}

    def close(self):
        self.ex.shutdown(wait=True)


class QueueWorkers:
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


def queue_worker(root, slots: int, poll: float = 0.5) -> None:
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


def make_workers(kind: str, n: int, queue_dir=None):
    """"local": n processes on this machine; "queue": the directory queue at queue_dir."""
    if kind == "local":
        return LocalWorkers(n)
    if kind == "queue":
        return QueueWorkers(queue_dir)
    raise ValueError(f"unknown workers '{kind}'")
