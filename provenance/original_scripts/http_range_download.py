#!/usr/bin/env python3
"""Resumable HTTP range downloader with byte-count and MD5 verification."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import threading
import time
import urllib.request
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--size", required=True, type=int)
    parser.add_argument("--md5", required=True)
    parser.add_argument("--workers", type=int, default=128)
    parser.add_argument("--chunk-mib", type=int, default=1)
    return parser.parse_args()


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        if args.output.stat().st_size == args.size and md5sum(args.output) == args.md5:
            print(f"Verified existing file: {args.output}")
            return
        raise SystemExit(f"Refusing to overwrite non-verified existing file: {args.output}")

    part = args.output.with_suffix(args.output.suffix + ".part")
    state_path = args.output.with_suffix(args.output.suffix + ".ranges.json")
    chunk_size = args.chunk_mib * 1024 * 1024
    ranges = [(start, min(start + chunk_size, args.size) - 1) for start in range(0, args.size, chunk_size)]
    completed: set[int] = set()
    if state_path.exists():
        state = json.loads(state_path.read_text())
        if state.get("url") != args.url or state.get("size") != args.size or state.get("chunk_size") != chunk_size:
            raise SystemExit(f"Resume state does not match requested download: {state_path}")
        completed = {int(value) for value in state.get("completed", [])}

    fd = os.open(part, os.O_RDWR | os.O_CREAT, 0o644)
    os.ftruncate(fd, args.size)
    lock = threading.Lock()

    def save_state() -> None:
        temporary = state_path.with_suffix(state_path.suffix + ".tmp")
        temporary.write_text(json.dumps({
            "url": args.url,
            "size": args.size,
            "chunk_size": chunk_size,
            "completed": sorted(completed),
        }))
        os.replace(temporary, state_path)

    def fetch(index: int) -> int:
        if index in completed:
            return index
        start, end = ranges[index]
        expected = end - start + 1
        last_error: Exception | None = None
        for attempt in range(10):
            try:
                request = urllib.request.Request(
                    args.url,
                    headers={
                        "Range": f"bytes={start}-{end}",
                        "User-Agent": "sleep-gut-project/PRJNA862187-metadata-verified",
                    },
                )
                with urllib.request.urlopen(request, timeout=180) as response:
                    content_range = response.headers.get("Content-Range", "")
                    payload = response.read()
                if len(payload) != expected:
                    raise IOError(f"range {index}: expected {expected} bytes, received {len(payload)}")
                if content_range and not content_range.startswith(f"bytes {start}-{end}/"):
                    raise IOError(f"range {index}: unexpected Content-Range {content_range}")
                os.pwrite(fd, payload, start)
                with lock:
                    completed.add(index)
                    if len(completed) % 25 == 0 or len(completed) == len(ranges):
                        save_state()
                        print(f"{args.output.name}: {len(completed)}/{len(ranges)} ranges", flush=True)
                return index
            except Exception as exc:
                last_error = exc
                time.sleep(min(30, 2 ** attempt))
        raise RuntimeError(f"failed range {index} after retries: {last_error}")

    try:
        pending = [index for index in range(len(ranges)) if index not in completed]
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for _ in pool.map(fetch, pending):
                pass
        os.fsync(fd)
    finally:
        os.close(fd)

    if part.stat().st_size != args.size:
        raise SystemExit(f"Size verification failed for {part}")
    observed = md5sum(part)
    if observed != args.md5:
        raise SystemExit(f"MD5 verification failed for {part}: expected {args.md5}, observed {observed}")
    os.replace(part, args.output)
    if state_path.exists():
        state_path.unlink()
    print(f"Completed and verified: {args.output}")


if __name__ == "__main__":
    main()
