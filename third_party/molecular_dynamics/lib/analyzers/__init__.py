"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
}
