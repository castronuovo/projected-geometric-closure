# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Vitantonio Castronuovo

"""Network-failure regression tests; no scientific calculations or downloads."""

import hashlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import HTTPError, URLError

from verified_download import download_verified


class InterruptedResponse(io.BytesIO):
    def read(self, size=-1):
        if self.tell():
            raise TimeoutError("interrupted transfer")
        return super().read(2)


class VerifiedDownloadTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "input.h5"
        self.data = b"verified input fixture"
        self.digest = hashlib.sha256(self.data).hexdigest()
        self.url = "https://example.invalid/input.h5"

    def assert_no_partial_files(self):
        self.assertEqual(list(self.path.parent.glob("*.part")), [])

    @patch("verified_download.time.sleep")
    @patch("verified_download.urlopen")
    def test_timeout_then_success(self, opener, sleep):
        opener.side_effect = [URLError(TimeoutError("connect")), io.BytesIO(self.data)]
        download_verified(self.url, self.path, self.digest, timeout=7)
        self.assertEqual(self.path.read_bytes(), self.data)
        self.assertEqual(opener.call_count, 2)
        opener.assert_called_with(self.url, timeout=7)
        sleep.assert_called_once_with(2)
        self.assert_no_partial_files()

    @patch("verified_download.time.sleep")
    @patch("verified_download.urlopen")
    def test_exhaustion_fails_without_cache_entry(self, opener, sleep):
        opener.side_effect = URLError(TimeoutError("connect"))
        with self.assertRaisesRegex(RuntimeError, "after 4 attempts"):
            download_verified(self.url, self.path, self.digest)
        self.assertEqual(opener.call_count, 4)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [2, 4, 8])
        self.assertFalse(self.path.exists())
        self.assert_no_partial_files()

    @patch("verified_download.time.sleep")
    @patch("verified_download.urlopen")
    def test_interrupted_transfer_is_removed_before_retry(self, opener, sleep):
        opener.side_effect = [InterruptedResponse(self.data), io.BytesIO(self.data)]
        download_verified(self.url, self.path, self.digest)
        self.assertEqual(self.path.read_bytes(), self.data)
        self.assertEqual(opener.call_count, 2)
        self.assert_no_partial_files()

    @patch("verified_download.urlopen")
    def test_wrong_hash_fails_without_retry_or_replacement(self, opener):
        self.path.write_bytes(b"old corrupt cache")
        opener.return_value = io.BytesIO(b"unexpected server content")
        with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
            download_verified(self.url, self.path, self.digest)
        self.assertEqual(opener.call_count, 1)
        self.assertEqual(self.path.read_bytes(), b"old corrupt cache")
        self.assert_no_partial_files()

    @patch("verified_download.urlopen")
    def test_valid_cache_avoids_network(self, opener):
        self.path.write_bytes(self.data)
        download_verified(self.url, self.path, self.digest)
        opener.assert_not_called()

    @patch("verified_download.urlopen")
    def test_corrupt_legacy_cache_is_replaced_only_after_verification(self, opener):
        self.path.write_bytes(b"partial old download")
        opener.return_value = io.BytesIO(self.data)
        download_verified(self.url, self.path, self.digest)
        self.assertEqual(self.path.read_bytes(), self.data)
        self.assert_no_partial_files()

    @patch("verified_download.urlopen")
    def test_permanent_http_error_is_not_retried(self, opener):
        body = io.BytesIO(b"missing")
        opener.side_effect = HTTPError(self.url, 404, "missing", {}, body)
        with self.assertRaises(HTTPError):
            download_verified(self.url, self.path, self.digest)
        self.assertEqual(opener.call_count, 1)
        self.assertTrue(body.closed)
        self.assert_no_partial_files()

    @patch("verified_download.time.sleep")
    @patch("verified_download.urlopen")
    def test_transient_http_error_is_retried(self, opener, sleep):
        opener.side_effect = [HTTPError(self.url, 503, "unavailable", {}, None),
                              io.BytesIO(self.data)]
        download_verified(self.url, self.path, self.digest)
        self.assertEqual(self.path.read_bytes(), self.data)
        self.assertEqual(opener.call_count, 2)


if __name__ == "__main__":
    unittest.main()
