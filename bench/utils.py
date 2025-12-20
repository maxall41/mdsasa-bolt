from pathlib import Path

from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis

from mdsasa_bolt.analysis import SASAAnalysis as SASAAnalysisBolt


def save(analysis: SASAAnalysis | SASAAnalysisBolt, filename: str) -> None:
    """Save analysis results by dumping frames as PDB with SASA in B-factor column."""
    with Path.open("bench/" / Path(filename), "w") as f:
        f.writelines(f"{frame}\n" for frame in analysis.results.total_area)
