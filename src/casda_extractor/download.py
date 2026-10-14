"""Parallel, retry-aware CASDA file downloader (shared by lightcurves and images)."""

from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)


def _fix_url(url: str) -> str:
    """Escape characters that break S3 / Squid proxy requests."""
    return url.replace("+", "%2B").replace(" ", "%20").replace('"', "%22")


def _download_single(
    url: str,
    savedir: str,
    casda_obj,
    file_num: int,
    total: int,
    max_retries: int = 3,
    retry_delay: float = 5.0,
) -> tuple[bool, str, str | None]:
    """
    Download one staged URL, retrying on gateway / timeout errors.

    Returns
    -------
    (success, filename_or_url, error_message)
    """
    escaped = _fix_url(url)
    for attempt in range(1, max_retries + 1):
        try:
            filelist = casda_obj.download_files([escaped], savedir=savedir)
            filename = filelist[0] if filelist else "unknown"
            print(f"  [{file_num}/{total}] Downloaded: {filename.split('/')[-1]}")
            return True, filename, None
        except Exception as exc:  # noqa: BLE001
            msg = str(exc)
            retryable = any(k in msg.lower() for k in ("timed out", "504", "502", "503"))
            if retryable and attempt < max_retries:
                snippet = msg[:60]
                print(
                    f"  [{file_num}/{total}] Attempt {attempt} "
                    f"failed ({snippet}…), retry in {retry_delay}s"
                )
                time.sleep(retry_delay)
            else:
                print(f"  [{file_num}/{total}] ERROR (attempt {attempt}/{max_retries}): {exc}")
                return False, escaped, msg
    return False, escaped, "Max retries exceeded"


def parallel_download(
    url_list: list[str],
    savedir: str,
    casda_obj,
    max_workers: int = 4,
    serial_retry: bool = True,
) -> list[str]:
    """
    Download multiple staged URLs in parallel, with an optional serial retry pass.

    Parameters
    ----------
    url_list : list[str]
        Staged download URLs from ``casda.stage_data()``.
    savedir : str
        Local directory to save files into.
    casda_obj : Casda
        Authenticated Casda instance.
    max_workers : int
        Parallel download threads.
    serial_retry : bool
        Retry failed downloads serially after the parallel pass.

    Returns
    -------
    list[str]
        Paths of successfully downloaded files.
    """
    downloaded: list[str] = []
    failed: list[tuple[str, str]] = []
    total = len(url_list)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_download_single, url, savedir, casda_obj, i + 1, total): url
            for i, url in enumerate(url_list)
        }
        for fut in as_completed(futures):
            try:
                ok, name, err = fut.result()
                if ok:
                    downloaded.append(name)
                else:
                    failed.append((futures[fut], err or ""))
            except Exception as exc:  # noqa: BLE001
                failed.append((futures[fut], str(exc)))

    print(f"\nDownload summary: {len(downloaded)}/{total} successful")

    if serial_retry and failed:
        print(f"\nRetrying {len(failed)} failed download(s) serially…")
        still_failed = []
        for url, _ in failed:
            ok, name, err = _download_single(url, savedir, casda_obj, 1, 1, max_retries=3, retry_delay=10)
            if ok:
                downloaded.append(name)
            else:
                still_failed.append((url, err or ""))
        failed = still_failed
        print(f"After retry: {len(downloaded)}/{total} successful")

    if failed:
        print(f"Failed ({len(failed)}):")
        for url, err in failed:
            print(f"  {url}: {err}")

    return downloaded


def stage_in_chunks(
    files_to_stage: list,
    all_data,
    casda_obj,
    chunk_size: int = 50,
    filename_col: str = "filename",
) -> list[str]:
    """
    Stage CASDA files in chunks (avoids URL-length limits), returning a deduplicated URL list.

    Parameters
    ----------
    files_to_stage : list
        Filenames to stage.
    all_data : pd.DataFrame
        Metadata table containing ``filename_col``.
    casda_obj : Casda
        Authenticated Casda instance.
    chunk_size : int
        Files to stage per batch.
    filename_col : str
        Column name for filenames.

    Returns
    -------
    list[str]
        Staged download URLs (deduplicated, checksums excluded).
    """
    url_list: list[str] = []
    n_chunks = (len(files_to_stage) + chunk_size - 1) // chunk_size

    for i in range(0, len(files_to_stage), chunk_size):
        chunk = files_to_stage[i : i + chunk_size]
        print(f"  Staging chunk {i // chunk_size + 1}/{n_chunks} ({len(chunk)} files)…")
        chunk_data = all_data[all_data[filename_col].isin(chunk)]
        try:
            chunk_table = chunk_data.to_records(index=False)
            urls = casda_obj.stage_data(chunk_table)
            for u in urls:
                if "checksum" not in u and u not in url_list:
                    url_list.append(u)
        except Exception as exc:  # noqa: BLE001
            print(f"  Chunk staging failed: {exc}. Falling back to individual staging…")
            for fname in chunk:
                try:
                    fd = all_data[all_data[filename_col] == fname]
                    ft = fd.to_records(index=False)
                    for u in casda_obj.stage_data(ft):
                        if "checksum" not in u and u not in url_list:
                            url_list.append(u)
                except Exception as exc2:  # noqa: BLE001
                    print(f"  Failed to stage {fname}: {exc2}")

    return url_list
