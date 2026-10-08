"""Rebuild and run both E6/E7/E8 terminal algorithms, then cross-check outputs.

Run from any directory: python3 research/en_e8/reproduce.py
Requires a C++17 compiler and Python standard library only. The final
cross-check also reads the independent certificates already in this project.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix="e8-terminal-reproduce-") as tmp:
        for stem, mode in [("parabolic_terminals", "flat"), ("recursive_terminals", "recursive")]:
            exe = Path(tmp) / stem
            subprocess.run(["c++", "-std=c++17", "-O3", str(HERE / (stem + ".cpp")), "-o", str(exe)], check=True)
            for n in (6, 7, 8):
                name = f"e{n}-recursive" if mode == "recursive" else ("e8-terminals" if n == 8 else f"e{n}-validation")
                result = subprocess.run([str(exe), str(n)], capture_output=True, text=True, check=True)
                data = json.loads(result.stdout)
                assert data["complete"]
                (HERE / (name + ".json")).write_text(result.stdout)
                (HERE / (name + ".log")).write_text(result.stderr)
        result = subprocess.run(["python3", str(HERE / "verify_outputs.py")], capture_output=True, text=True, check=True)
        (HERE / "cross-checks.json").write_text(result.stdout)
        print(result.stdout, end="")
    sources = [HERE / name for name in ["parabolic_terminals.cpp", "recursive_terminals.cpp", "verify_outputs.py", "reproduce.py"]]
    manifest = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (HERE / "source-hashes.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
