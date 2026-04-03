"""Download PDB files from RCSB for validation."""
import os
import sys
import urllib.request
import concurrent.futures
import time

def download_pdb(pdb_id, out_dir):
    """Download a single PDB file. Returns (pdb_id, success, error)."""
    out_path = os.path.join(out_dir, f"{pdb_id.lower()}.pdb")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return (pdb_id, True, "exists")

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        urllib.request.urlretrieve(url, out_path)
        return (pdb_id, True, None)
    except Exception as e:
        # Try mmCIF as fallback (some entries don't have PDB format)
        cif_path = os.path.join(out_dir, f"{pdb_id.lower()}.cif")
        url_cif = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        try:
            urllib.request.urlretrieve(url_cif, cif_path)
            return (pdb_id, True, "cif_fallback")
        except Exception as e2:
            return (pdb_id, False, str(e2)[:100])

def main():
    id_file = sys.argv[1] if len(sys.argv) > 1 else "validation/pdb_ids_10k.txt"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "validation/pdbs_10k/"
    max_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    os.makedirs(out_dir, exist_ok=True)

    with open(id_file) as f:
        pdb_ids = [line.strip() for line in f if line.strip()]

    print(f"Downloading {len(pdb_ids)} PDB files to {out_dir} ({max_workers} workers)")

    done = 0
    failed = 0
    skipped = 0
    start = time.time()

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(download_pdb, pid, out_dir): pid for pid in pdb_ids}
        for future in concurrent.futures.as_completed(futures):
            pdb_id, success, error = future.result()
            done += 1
            if not success:
                failed += 1
            elif error == "exists":
                skipped += 1

            if done % 500 == 0 or done == len(pdb_ids):
                elapsed = time.time() - start
                rate = done / elapsed
                eta = (len(pdb_ids) - done) / rate if rate > 0 else 0
                print(f"  [{done}/{len(pdb_ids)}] {failed} failed, {skipped} cached, "
                      f"{rate:.0f}/s, ETA {eta:.0f}s")

    elapsed = time.time() - start
    print(f"\nDone in {elapsed:.0f}s: {done - failed} downloaded, {failed} failed")

if __name__ == "__main__":
    main()
