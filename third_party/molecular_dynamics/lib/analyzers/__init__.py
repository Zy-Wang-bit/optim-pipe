"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer
from .rmsf import RMSFAnalyzer
from .hbond import HBondAnalyzer
from .sasa import SASAAnalyzer
from .contacts import ContactsAnalyzer
from .mmpbsa import MMPBSAAnalyzer
from .pairwise_dist import PairwiseDistAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
    "rmsf": RMSFAnalyzer,
    "hbond": HBondAnalyzer,
    "sasa": SASAAnalyzer,
    "contacts": ContactsAnalyzer,
    "mmpbsa": MMPBSAAnalyzer,
    "pairwise_dist": PairwiseDistAnalyzer,
}
