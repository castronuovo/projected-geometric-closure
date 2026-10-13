# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Download immutable inputs with bounded retries and verified atomic writes."""

import hashlib
from http.client import IncompleteRead
import os
from pathlib import Path
import shutil
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_verified(url: str, path: Path, expected_hash: str,
                      *, timeout: float = 45.0, attempts: int = 4) -> Path:
    """Reuse only valid cache entries; publish a download only after hashing.

    ``timeout`` limits blocking socket operations, not total elapsed transfer
    time. Transient network failures and HTTP 408/429/selected 5xx are retried with
    bounded backoff. A digest mismatch or permanent HTTP error fails directly.
    """
    if timeout <= 0 or attempts < 1:
        raise ValueError("timeout must be positive and attempts at least one")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and sha256(path) == expected_hash:
        return path

    for attempt in range(attempts):
        fd, temporary_name = tempfile.mkstemp(
            dir=path.parent, prefix=path.name + ".", suffix=".part"
        )
        temporary = Path(temporary_name)
        try:
            with os.fdopen(fd, "wb") as output:
                with urlopen(url, timeout=timeout) as response:
                    shutil.copyfileobj(response, output)
            observed_hash = sha256(temporary)
            if observed_hash != expected_hash:
                raise RuntimeError(
                    f"SHA-256 mismatch for {path.name}: {observed_hash}; "
                    f"expected {expected_hash}"
                )
            temporary.replace(path)
            return path
        except HTTPError as error:
            if error.fp is not None:
                error.close()
            if error.code not in (408, 429, 500, 502, 503, 504):
                raise
            last_error = error
        except (URLError, TimeoutError, ConnectionError, IncompleteRead) as error:
            last_error = error
        finally:
            temporary.unlink(missing_ok=True)

        if attempt + 1 < attempts:
            delay = 2 ** (attempt + 1)
            print(f"Download attempt {attempt + 1}/{attempts} failed for "
                  f"{path.name}: {last_error}; retrying in {delay}s", flush=True)
            time.sleep(delay)

    raise RuntimeError(
        f"Download failed for {path.name} after {attempts} attempts "
        f"(socket timeout {timeout}s): {url}"
    ) from last_error
