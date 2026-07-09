# Make the emperor package importable in tests regardless of install/venv state.
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))
