# MD 分析模块实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an independent MD module that automates GROMACS execution (with CpHMD support) and provides a comprehensive trajectory analysis toolkit for evaluating structural stability, interface integrity, pH sensitivity, and binding affinity.

**Architecture:** Python wrapper around GROMACS commands with Jinja2-templated MDP files for flexible simulation parameters. Analysis uses MDAnalysis for trajectory processing. All config in `md_config.yaml`, output directory configurable (default: `experiments/<system>/<round>/md/`).

**Tech Stack:** Python 3.11, GROMACS (≥2022), Jinja2 (installed), MDAnalysis (to install), gmx_MMPBSA (optional, to install), PyYAML, pandas, numpy

**Spec:** `docs/superpowers/specs/2026-04-07-md-analysis-module-design.md`

---

## File Structure

### New files to create:

```
third_party/molecular_dynamics/
├── configs/
│   └── md_config.yaml              # MD module configuration
├── mdp_templates/
│   ├── minimization.mdp.j2         # Energy minimization template
│   ├── nvt.mdp.j2                  # NVT equilibration template
│   ├── npt.mdp.j2                  # NPT equilibration template
│   └── production.mdp.j2           # Production MD template (standard + CpHMD)
├── run_md.py                        # Execution entry point
├── analyze_trajectory.py            # Analysis entry point
├── compare_ph.py                    # Dual-pH comparison
├── lib/
│   ├── __init__.py
│   ├── config.py                    # Config loading + validation
│   ├── gromacs_wrapper.py           # GROMACS command wrapper
│   ├── protonation.py               # His protonation state management
│   ├── analyzers/
│   │   ├── __init__.py
│   │   ├── base.py                  # Base analyzer class
│   │   ├── rmsd.py                  # RMSD time series + convergence
│   │   ├── rmsf.py                  # Per-residue flexibility
│   │   ├── hbond.py                 # Interface H-bond occupancy
│   │   ├── sasa.py                  # Buried SASA
│   │   ├── contacts.py              # Interface contacts
│   │   └── mmpbsa.py                # MM-GBSA (optional)
│   └── report.py                    # Summary report generation
├── tests/
│   ├── __init__.py
│   ├── test_config.py
│   ├── test_gromacs_wrapper.py
│   ├── test_protonation.py
│   ├── test_analyzers.py
│   └── test_report.py
└── legacy/                          # Move existing scripts here
    ├── autoGromacs.sh
    ├── true-autoGromacs.sh
    └── show_xvg.py
```

### Existing files to move:

- `third_party/molecular_dynamics/autoGromacs.sh` → `legacy/`
- `third_party/molecular_dynamics/true-autoGromacs.sh` → `legacy/`
- `third_party/molecular_dynamics/show_xvg.py` → `legacy/`
- `third_party/molecular_dynamics/1.minimization.mdp` → remove after templates created
- `third_party/molecular_dynamics/2.1equilibration_NVT.mdp` → remove after templates created
- `third_party/molecular_dynamics/2.2equilibration_NPT.mdp` → remove after templates created
- `third_party/molecular_dynamics/3.production.mdp` → remove after templates created

---

## Task 1: Project Setup — Directory Structure, Config, Dependencies

**Files:**
- Create: `third_party/molecular_dynamics/configs/md_config.yaml`
- Create: `third_party/molecular_dynamics/lib/__init__.py`
- Create: `third_party/molecular_dynamics/lib/config.py`
- Create: `third_party/molecular_dynamics/lib/analyzers/__init__.py`
- Create: `third_party/molecular_dynamics/tests/__init__.py`
- Create: `third_party/molecular_dynamics/tests/test_config.py`
- Move: existing scripts → `legacy/`

- [ ] **Step 1: Create directory structure and move legacy files**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
mkdir -p configs mdp_templates lib/analyzers tests legacy

# Move existing scripts to legacy/
mv autoGromacs.sh true-autoGromacs.sh show_xvg.py legacy/

# Move old MDP files (will be replaced by templates)
mv 1.minimization.mdp 2.1equilibration_NVT.mdp 2.2equilibration_NPT.mdp 3.production.mdp legacy/

# Create __init__.py files
touch lib/__init__.py lib/analyzers/__init__.py tests/__init__.py
```

- [ ] **Step 2: Write md_config.yaml**

```yaml
# third_party/molecular_dynamics/configs/md_config.yaml
# MD 模块配置文件

md:
  # GROMACS 执行参数
  gmx_executable: "gmx"           # gmx 可执行文件路径，可改为绝对路径
  force_field: "amber99sb-ildn"    # 力场名称（传给 -ff 参数）
  water_model: "tip3p"             # 水模型名称（传给 -water 参数）
  box_type: "cubic"
  box_distance: 1.0               # nm, 蛋白到盒子边缘的最小距离
  ion_concentration: 0.15          # mol/L
  temperature: 310                 # K (37°C)
  max_warnings: 2

  stages:
    minimization:
      max_steps: 5000
      emtol: 1000.0                # kJ/mol/nm
    nvt:
      duration_ps: 200
      dt_ps: 0.001
    npt:
      duration_ps: 1000
      dt_ps: 0.001
    production:
      duration_ns: 50
      dt_ps: 0.002
      save_interval_ps: 100        # 轨迹保存间隔

  cphmd:
    enabled: false                  # CpHMD 默认关闭（需 GROMACS ≥2022）
    ph_values: [7.4, 6.0]
    titratable_residues: []         # 空 = 自动检测所有 His
    lambda_dynamics:
      calibration_ph: 7.0
      barrier_height: 5.0           # kJ/mol
      update_nst: 100               # λ 更新间隔 (steps)

  gpu:
    device_id: 0
    ntmpi: 1
    ntomp: 5

  output_dir: "experiments/1E62/R3/md"

analysis:
  antibody_chains: ["A", "B"]
  antigen_chains: ["C"]

  cdr_regions:
    H1: { chain: "A", start: 26, end: 35 }
    H2: { chain: "A", start: 50, end: 65 }
    H3: { chain: "A", start: 95, end: 102 }
    L1: { chain: "B", start: 24, end: 34 }
    L2: { chain: "B", start: 50, end: 56 }
    L3: { chain: "B", start: 89, end: 97 }

  convergence:
    window_ns: 5.0
    variance_threshold: 0.05       # Å²

  contact_distance: 4.5             # Å

  mmpbsa:
    enabled: false
    method: "GBSA"
    frame_interval_ps: 100
```

- [ ] **Step 3: Write config.py — config loading and validation**

```python
# third_party/molecular_dynamics/lib/config.py
"""MD 模块配置加载与校验。"""

from pathlib import Path
from typing import Any

import yaml

_MODULE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = _MODULE_ROOT / "configs" / "md_config.yaml"


def load_config(path: Path | str | None = None) -> dict[str, Any]:
    """加载并校验 md_config.yaml。

    Parameters
    ----------
    path : 配置文件路径，默认 configs/md_config.yaml

    Returns
    -------
    dict  完整配置字典
    """
    path = Path(path) if path else DEFAULT_CONFIG
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    _validate(cfg)
    return cfg


def _validate(cfg: dict) -> None:
    """基本校验：必需字段存在且类型正确。"""
    md = cfg.get("md")
    if md is None:
        raise ValueError("Config missing top-level 'md' section")

    required_md = ["force_field", "water_model", "temperature", "stages"]
    for key in required_md:
        if key not in md:
            raise ValueError(f"Config md.{key} is required")

    stages = md["stages"]
    for stage in ("minimization", "nvt", "npt", "production"):
        if stage not in stages:
            raise ValueError(f"Config md.stages.{stage} is required")

    analysis = cfg.get("analysis")
    if analysis is None:
        raise ValueError("Config missing top-level 'analysis' section")

    for key in ("antibody_chains", "antigen_chains", "cdr_regions"):
        if key not in analysis:
            raise ValueError(f"Config analysis.{key} is required")


def get_production_nsteps(cfg: dict) -> int:
    """从配置计算 production 阶段的总步数。"""
    prod = cfg["md"]["stages"]["production"]
    duration_ps = prod["duration_ns"] * 1000
    dt_ps = prod["dt_ps"]
    return int(duration_ps / dt_ps)


def get_save_nsteps(cfg: dict) -> int:
    """从配置计算轨迹保存间隔的步数。"""
    prod = cfg["md"]["stages"]["production"]
    return int(prod["save_interval_ps"] / prod["dt_ps"])
```

- [ ] **Step 4: Write test for config loading**

```python
# third_party/molecular_dynamics/tests/test_config.py
"""Tests for config loading and validation."""

import tempfile
from pathlib import Path

import pytest
import yaml

from lib.config import load_config, get_production_nsteps, get_save_nsteps


def test_load_default_config():
    cfg = load_config()
    assert cfg["md"]["force_field"] == "amber99sb-ildn"
    assert cfg["md"]["temperature"] == 310
    assert "H1" in cfg["analysis"]["cdr_regions"]


def test_load_missing_file():
    with pytest.raises(FileNotFoundError):
        load_config("/nonexistent/path.yaml")


def test_validation_missing_md_section():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump({"analysis": {"antibody_chains": ["A"]}}, f)
        f.flush()
        with pytest.raises(ValueError, match="'md' section"):
            load_config(f.name)


def test_get_production_nsteps():
    cfg = load_config()
    nsteps = get_production_nsteps(cfg)
    # 50 ns / 0.002 ps = 25,000,000
    assert nsteps == 25_000_000


def test_get_save_nsteps():
    cfg = load_config()
    nsteps = get_save_nsteps(cfg)
    # 100 ps / 0.002 ps = 50,000
    assert nsteps == 50_000
```

- [ ] **Step 5: Run tests**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python -m pytest tests/test_config.py -v
```

Expected: All 5 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add third_party/molecular_dynamics/
git commit -m "feat(md): scaffold MD module with config, move legacy scripts"
```

---

## Task 2: MDP Jinja2 Templates

**Files:**
- Create: `third_party/molecular_dynamics/mdp_templates/minimization.mdp.j2`
- Create: `third_party/molecular_dynamics/mdp_templates/nvt.mdp.j2`
- Create: `third_party/molecular_dynamics/mdp_templates/npt.mdp.j2`
- Create: `third_party/molecular_dynamics/mdp_templates/production.mdp.j2`

Templates are parameterized versions of the existing MDP files. Variables come from `md_config.yaml`.

- [ ] **Step 1: Write minimization.mdp.j2**

```jinja2
; third_party/molecular_dynamics/mdp_templates/minimization.mdp.j2
; Energy minimization — generated from template
define                  = -DPOSRES
integrator              = steep
emtol                   = {{ emtol }}
nsteps                  = {{ max_steps }}
nstlist                 = 10
cutoff-scheme           = Verlet
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = pme
rcoulomb                = 1.2
;
constraints             = h-bonds
constraint_algorithm    = LINCS
```

- [ ] **Step 2: Write nvt.mdp.j2**

```jinja2
; third_party/molecular_dynamics/mdp_templates/nvt.mdp.j2
; NVT equilibration — generated from template
define                  = -DPOSRES
integrator              = md
dt                      = {{ dt }}
nsteps                  = {{ nsteps }}
nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = pme
rcoulomb                = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = Protein Non-Protein
tau_t                   = 1.0 1.0
ref_t                   = {{ temperature }} {{ temperature }}
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
gen-vel                 = yes
gen-temp                = {{ temperature }}
gen-seed                = -1
```

- [ ] **Step 3: Write npt.mdp.j2**

```jinja2
; third_party/molecular_dynamics/mdp_templates/npt.mdp.j2
; NPT equilibration — generated from template
define                  = -DPOSRES
integrator              = md
dt                      = {{ dt }}
nsteps                  = {{ nsteps }}
nstxtcout               = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = pme
rcoulomb                = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = Protein Non-Protein
tau_t                   = 1.0 1.0
ref_t                   = {{ temperature }} {{ temperature }}
;
Pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = 1.0
compressibility         = 4.5e-5
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
gen-vel                 = yes
gen-temp                = {{ temperature }}
gen-seed                = -1
```

- [ ] **Step 4: Write production.mdp.j2**

```jinja2
; third_party/molecular_dynamics/mdp_templates/production.mdp.j2
; Production MD — generated from template
integrator              = md
dt                      = {{ dt }}
nsteps                  = {{ nsteps }}
;
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 50000
nstenergy               = 10000
nstxout-compressed      = {{ save_nsteps }}
compressed-x-grps       = System
;
cutoff-scheme           = Verlet
nstlist                 = 20
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
rcoulomb                = 1.2
coulombtype             = pme
;
tcoupl                  = Nose-Hoover
tc_grps                 = Protein Non-Protein
tau_t                   = 1.0 1.0
ref_t                   = {{ temperature }} {{ temperature }}
;
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = System
{% if cphmd_enabled %}
;
; --- Constant pH MD (λ-dynamics) ---
lambda-dynamics          = yes
lambda-dynamics-simulation-ph = {{ cphmd_ph }}
lambda-dynamics-lambda-group-count = {{ cphmd_n_groups }}
lambda-dynamics-update-nst = {{ cphmd_update_nst }}
{% for group in cphmd_groups %}
lambda-dynamics-group-type{{ loop.index }} = {{ group.type }}
lambda-dynamics-atom-set{{ loop.index }}-name = {{ group.name }}
lambda-dynamics-atom-set{{ loop.index }}-index-group-name = {{ group.index_group }}
lambda-dynamics-atom-set{{ loop.index }}-barrier = {{ group.barrier }}
lambda-dynamics-atom-set{{ loop.index }}-init-lambda = {{ group.init_lambda }}
lambda-dynamics-charge-constraints = {{ group.charge_constraints | default('no') }}
{% endfor %}
{% endif %}
```

- [ ] **Step 5: Commit**

```bash
git add third_party/molecular_dynamics/mdp_templates/
git commit -m "feat(md): add Jinja2 MDP templates for all simulation stages"
```

---

## Task 3: GROMACS Wrapper — Subprocess Management

**Files:**
- Create: `third_party/molecular_dynamics/lib/gromacs_wrapper.py`
- Create: `third_party/molecular_dynamics/tests/test_gromacs_wrapper.py`

This is the core automation layer. It wraps each GROMACS command as a Python function, uses `-ff`/`-water` name-based selection instead of interactive indices, and renders MDP templates.

- [ ] **Step 1: Write gromacs_wrapper.py**

```python
# third_party/molecular_dynamics/lib/gromacs_wrapper.py
"""GROMACS 命令封装：解决 pdb2gmx 交互选择问题，实现全自动 MD 流程。"""

import logging
import subprocess
from pathlib import Path
from typing import Any

import jinja2

from .config import get_production_nsteps, get_save_nsteps

logger = logging.getLogger(__name__)

_TEMPLATE_DIR = Path(__file__).resolve().parent.parent / "mdp_templates"


class GromacsError(RuntimeError):
    """GROMACS 命令执行失败。"""

    def __init__(self, cmd: str, returncode: int, stderr: str):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(f"gmx {cmd} failed (rc={returncode}):\n{stderr[-2000:]}")


class GromacsWrapper:
    """GROMACS 命令的 Python 封装。

    Parameters
    ----------
    cfg : dict
        从 md_config.yaml 加载的完整配置。
    work_dir : Path
        工作目录（所有 GROMACS 命令在此目录下执行）。
    """

    def __init__(self, cfg: dict, work_dir: Path):
        self.cfg = cfg
        self.md = cfg["md"]
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.gmx = self.md.get("gmx_executable", "gmx")
        self._jinja_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(str(_TEMPLATE_DIR)),
            undefined=jinja2.StrictUndefined,
        )

    # ---- 底层执行 ----

    def _run(
        self,
        args: list[str],
        stdin_text: str | None = None,
        check: bool = True,
    ) -> subprocess.CompletedProcess:
        """执行 GROMACS 命令。"""
        cmd = [self.gmx] + args
        logger.info("Running: %s", " ".join(cmd))
        result = subprocess.run(
            cmd,
            cwd=self.work_dir,
            input=stdin_text,
            capture_output=True,
            text=True,
            timeout=3600 * 24,  # 24h max per command
        )
        if check and result.returncode != 0:
            raise GromacsError(args[0], result.returncode, result.stderr)
        return result

    def _render_mdp(
        self, template_name: str, output_name: str, **kwargs: Any
    ) -> Path:
        """渲染 Jinja2 MDP 模板并写入工作目录。"""
        tmpl = self._jinja_env.get_template(template_name)
        content = tmpl.render(**kwargs)
        out = self.work_dir / output_name
        out.write_text(content)
        logger.info("Rendered MDP: %s", out)
        return out

    # ---- MD 流程步骤 ----

    def pdb2gmx(self, pdb_path: Path) -> Path:
        """拓扑生成。使用 -ff/-water 名称指定，不依赖交互索引。

        Returns processed .gro file path.
        """
        out_gro = self.work_dir / "processed.gro"
        args = [
            "pdb2gmx",
            "-f", str(pdb_path),
            "-o", str(out_gro),
            "-p", str(self.work_dir / "topol.top"),
            "-ff", self.md["force_field"],
            "-water", self.md["water_model"],
            "-ignh",
        ]
        self._run(args)
        return out_gro

    def editconf(self, gro_path: Path) -> Path:
        """创建模拟盒子。"""
        out_gro = self.work_dir / "box.gro"
        self._run([
            "editconf",
            "-f", str(gro_path),
            "-o", str(out_gro),
            "-d", str(self.md["box_distance"]),
            "-bt", self.md["box_type"],
        ])
        return out_gro

    def solvate(self, gro_path: Path) -> Path:
        """填充溶剂。"""
        out_gro = self.work_dir / "solvated.gro"
        self._run([
            "solvate",
            "-cp", str(gro_path),
            "-o", str(out_gro),
            "-p", str(self.work_dir / "topol.top"),
        ])
        return out_gro

    def add_ions(self, gro_path: Path) -> Path:
        """添加离子至中性 + 指定浓度。"""
        # 先 grompp 生成 tpr
        minim_mdp = self._render_minimization_mdp()
        tpr = self.work_dir / "ions.tpr"
        self._run([
            "grompp",
            "-f", str(minim_mdp),
            "-c", str(gro_path),
            "-r", str(gro_path),
            "-p", str(self.work_dir / "topol.top"),
            "-o", str(tpr),
            "-maxwarn", str(self.md["max_warnings"]),
        ])
        # genion: 选择 SOL 组来替换水分子
        out_gro = self.work_dir / "ions.gro"
        self._run(
            [
                "genion",
                "-s", str(tpr),
                "-p", str(self.work_dir / "topol.top"),
                "-o", str(out_gro),
                "-neutral",
                "-conc", str(self.md["ion_concentration"]),
            ],
            stdin_text="SOL\n",
        )
        return out_gro

    def minimization(self, gro_path: Path) -> Path:
        """能量最小化。"""
        mdp = self._render_minimization_mdp()
        tpr = self.work_dir / "em.tpr"
        self._run([
            "grompp",
            "-f", str(mdp),
            "-c", str(gro_path),
            "-r", str(gro_path),
            "-p", str(self.work_dir / "topol.top"),
            "-o", str(tpr),
            "-maxwarn", str(self.md["max_warnings"]),
        ])
        self._run([
            "mdrun",
            "-v", "-deffnm", "em",
            "-ntomp", str(self.md["gpu"]["ntomp"]),
            "-ntmpi", str(self.md["gpu"]["ntmpi"]),
            "-gpu_id", str(self.md["gpu"]["device_id"]),
        ])
        return self.work_dir / "em.gro"

    def nvt_equilibration(self, gro_path: Path) -> Path:
        """NVT 平衡。"""
        nvt_cfg = self.md["stages"]["nvt"]
        nsteps = int(nvt_cfg["duration_ps"] / nvt_cfg["dt_ps"])
        mdp = self._render_mdp(
            "nvt.mdp.j2", "nvt.mdp",
            dt=nvt_cfg["dt_ps"],
            nsteps=nsteps,
            temperature=self.md["temperature"],
        )
        tpr = self.work_dir / "NVT.tpr"
        self._run([
            "grompp",
            "-f", str(mdp),
            "-c", str(gro_path),
            "-r", str(gro_path),
            "-p", str(self.work_dir / "topol.top"),
            "-o", str(tpr),
            "-maxwarn", str(self.md["max_warnings"]),
        ])
        self._run([
            "mdrun", "-v", "-deffnm", "NVT",
            "-ntomp", str(self.md["gpu"]["ntomp"]),
            "-ntmpi", str(self.md["gpu"]["ntmpi"]),
            "-gpu_id", str(self.md["gpu"]["device_id"]),
        ])
        return self.work_dir / "NVT.gro"

    def npt_equilibration(self, gro_path: Path) -> Path:
        """NPT 平衡。"""
        npt_cfg = self.md["stages"]["npt"]
        nsteps = int(npt_cfg["duration_ps"] / npt_cfg["dt_ps"])
        mdp = self._render_mdp(
            "npt.mdp.j2", "npt.mdp",
            dt=npt_cfg["dt_ps"],
            nsteps=nsteps,
            temperature=self.md["temperature"],
        )
        tpr = self.work_dir / "NPT.tpr"
        self._run([
            "grompp",
            "-f", str(mdp),
            "-c", str(gro_path),
            "-r", str(gro_path),
            "-p", str(self.work_dir / "topol.top"),
            "-o", str(tpr),
            "-maxwarn", str(self.md["max_warnings"]),
        ])
        self._run([
            "mdrun", "-v", "-deffnm", "NPT",
            "-ntomp", str(self.md["gpu"]["ntomp"]),
            "-ntmpi", str(self.md["gpu"]["ntmpi"]),
            "-gpu_id", str(self.md["gpu"]["device_id"]),
        ])
        return self.work_dir / "NPT.gro"

    def production(self, gro_path: Path, ph: float | None = None) -> Path:
        """Production MD。

        Parameters
        ----------
        gro_path : NPT 平衡后的 .gro 文件
        ph : CpHMD 的目标 pH。None 表示标准 MD。

        Returns
        -------
        Path  production .tpr 文件路径
        """
        prod_cfg = self.md["stages"]["production"]
        nsteps = get_production_nsteps(self.cfg)
        save_nsteps = get_save_nsteps(self.cfg)

        # CpHMD 参数
        cphmd_enabled = (
            ph is not None and self.md.get("cphmd", {}).get("enabled", False)
        )
        cphmd_kwargs = {}
        if cphmd_enabled:
            cphmd_cfg = self.md["cphmd"]
            cphmd_kwargs = {
                "cphmd_enabled": True,
                "cphmd_ph": ph,
                "cphmd_n_groups": 0,  # 由 protonation.py 设置
                "cphmd_update_nst": cphmd_cfg["lambda_dynamics"]["update_nst"],
                "cphmd_groups": [],   # 由 protonation.py 填充
            }

        mdp = self._render_mdp(
            "production.mdp.j2", "production.mdp",
            dt=prod_cfg["dt_ps"],
            nsteps=nsteps,
            save_nsteps=save_nsteps,
            temperature=self.md["temperature"],
            **cphmd_kwargs,
        )

        tpr = self.work_dir / "production.tpr"
        self._run([
            "grompp",
            "-f", str(mdp),
            "-c", str(gro_path),
            "-r", str(self.work_dir / "NVT.gro"),
            "-p", str(self.work_dir / "topol.top"),
            "-o", str(tpr),
            "-maxwarn", str(self.md["max_warnings"]),
        ])
        self._run([
            "mdrun",
            "-s", str(tpr),
            "-v",
            "-deffnm", "production",
            "-ntomp", str(self.md["gpu"]["ntomp"]),
            "-ntmpi", str(self.md["gpu"]["ntmpi"]),
            "-gpu_id", str(self.md["gpu"]["device_id"]),
        ])
        return tpr

    # ---- 完整流程 ----

    def run_full_pipeline(self, pdb_path: Path, ph: float | None = None) -> Path:
        """执行完整 MD 流程：PDB → production 轨迹。

        Returns production .tpr path.
        """
        logger.info("=== Starting MD pipeline for %s (pH=%s) ===", pdb_path.name, ph)
        gro = self.pdb2gmx(pdb_path)
        gro = self.editconf(gro)
        gro = self.solvate(gro)
        gro = self.add_ions(gro)
        gro = self.minimization(gro)
        gro = self.nvt_equilibration(gro)
        gro = self.npt_equilibration(gro)
        tpr = self.production(gro, ph=ph)
        logger.info("=== MD pipeline complete ===")
        return tpr

    # ---- 内部 helpers ----

    def _render_minimization_mdp(self) -> Path:
        """渲染 minimization MDP（被 add_ions 和 minimization 共用）。"""
        minim_cfg = self.md["stages"]["minimization"]
        return self._render_mdp(
            "minimization.mdp.j2", "minimization.mdp",
            max_steps=minim_cfg["max_steps"],
            emtol=minim_cfg["emtol"],
        )
```

- [ ] **Step 2: Write unit tests for wrapper (mock subprocess)**

```python
# third_party/molecular_dynamics/tests/test_gromacs_wrapper.py
"""Tests for GromacsWrapper — mock subprocess calls."""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from lib.config import load_config
from lib.gromacs_wrapper import GromacsWrapper, GromacsError


@pytest.fixture
def cfg():
    return load_config()


@pytest.fixture
def wrapper(cfg, tmp_path):
    return GromacsWrapper(cfg, tmp_path)


@patch("lib.gromacs_wrapper.subprocess.run")
def test_run_success(mock_run, wrapper):
    mock_run.return_value = subprocess.CompletedProcess(
        args=["gmx", "help"], returncode=0, stdout="ok", stderr=""
    )
    result = wrapper._run(["help"])
    assert result.returncode == 0
    mock_run.assert_called_once()


@patch("lib.gromacs_wrapper.subprocess.run")
def test_run_failure_raises(mock_run, wrapper):
    mock_run.return_value = subprocess.CompletedProcess(
        args=["gmx", "pdb2gmx"], returncode=1, stdout="", stderr="Fatal error"
    )
    with pytest.raises(GromacsError, match="pdb2gmx failed"):
        wrapper._run(["pdb2gmx", "-f", "test.pdb"])


def test_render_mdp_minimization(wrapper):
    mdp_path = wrapper._render_minimization_mdp()
    assert mdp_path.exists()
    content = mdp_path.read_text()
    assert "integrator" in content
    assert "5000" in content  # max_steps from config


def test_render_mdp_production(wrapper):
    mdp_path = wrapper._render_mdp(
        "production.mdp.j2", "test_prod.mdp",
        dt=0.002, nsteps=25000000, save_nsteps=50000,
        temperature=310, cphmd_enabled=False,
    )
    content = mdp_path.read_text()
    assert "25000000" in content
    assert "lambda-dynamics" not in content  # CpHMD disabled


def test_render_mdp_production_with_cphmd(wrapper):
    mdp_path = wrapper._render_mdp(
        "production.mdp.j2", "test_cphmd.mdp",
        dt=0.002, nsteps=25000000, save_nsteps=50000,
        temperature=310,
        cphmd_enabled=True, cphmd_ph=6.0, cphmd_n_groups=1,
        cphmd_update_nst=100,
        cphmd_groups=[{
            "type": "histidine",
            "name": "HIS_42",
            "index_group": "HIS_42",
            "barrier": 5.0,
            "init_lambda": 0.5,
        }],
    )
    content = mdp_path.read_text()
    assert "lambda-dynamics" in content
    assert "6.0" in content
```

- [ ] **Step 3: Run tests**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python -m pytest tests/test_gromacs_wrapper.py -v
```

Expected: All tests PASS.

- [ ] **Step 4: Commit**

```bash
git add third_party/molecular_dynamics/lib/gromacs_wrapper.py \
        third_party/molecular_dynamics/tests/test_gromacs_wrapper.py
git commit -m "feat(md): add GROMACS wrapper with name-based force field selection"
```

---

## Task 4: Protonation Manager (CpHMD Setup)

**Files:**
- Create: `third_party/molecular_dynamics/lib/protonation.py`
- Create: `third_party/molecular_dynamics/tests/test_protonation.py`

- [ ] **Step 1: Write protonation.py**

```python
# third_party/molecular_dynamics/lib/protonation.py
"""His 质子化态管理，用于 CpHMD λ-dynamics 配置。

职责：
1. 从 PDB 文件中检测所有 His 残基位置
2. 根据 config 过滤要参与 CpHMD 的残基列表
3. 生成 CpHMD 所需的 λ-dynamics 组参数
"""

import logging
from pathlib import Path

from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)


def detect_his_residues(pdb_path: Path) -> list[dict]:
    """从 PDB 文件中检测所有 His 残基。

    Returns
    -------
    list[dict]  每个元素 {"chain": str, "resid": int, "resname": str}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", str(pdb_path))
    his_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue  # skip HETATM
                if residue.resname in ("HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"):
                    his_list.append({
                        "chain": chain.id,
                        "resid": residue.id[1],
                        "resname": residue.resname,
                    })
        break  # only first model
    logger.info("Detected %d His residues in %s", len(his_list), pdb_path.name)
    return his_list


def filter_titratable(
    all_his: list[dict],
    config_residues: list[dict],
) -> list[dict]:
    """过滤参与 CpHMD 的 His 残基。

    Parameters
    ----------
    all_his : detect_his_residues 返回的列表
    config_residues : md_config.yaml 中 cphmd.titratable_residues
                     空列表 = 使用全部检测到的 His

    Returns
    -------
    list[dict]  过滤后的 His 列表
    """
    if not config_residues:
        return all_his
    keep = {(r["chain"], r["resid"]) for r in config_residues}
    return [h for h in all_his if (h["chain"], h["resid"]) in keep]


def build_cphmd_groups(
    titratable_his: list[dict],
    barrier_height: float = 5.0,
) -> list[dict]:
    """为每个可滴定 His 生成 CpHMD λ-dynamics 组参数。

    Parameters
    ----------
    titratable_his : filter_titratable 返回的列表
    barrier_height : λ 势垒高度 (kJ/mol)

    Returns
    -------
    list[dict]  CpHMD 组参数，可直接传给 production.mdp.j2 的 cphmd_groups
    """
    groups = []
    for his in titratable_his:
        name = f"HIS_{his['chain']}{his['resid']}"
        groups.append({
            "type": "histidine",
            "name": name,
            "index_group": name,
            "barrier": barrier_height,
            "init_lambda": 0.5,  # 从中间态开始
        })
    return groups
```

- [ ] **Step 2: Write tests**

```python
# third_party/molecular_dynamics/tests/test_protonation.py
"""Tests for protonation manager."""

import tempfile
from pathlib import Path

import pytest

from lib.protonation import detect_his_residues, filter_titratable, build_cphmd_groups


# Minimal PDB with 2 His residues for testing
_MINI_PDB = """\
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   1.000   1.000  1.00  0.00           C
ATOM      3  N   HIS A  42       3.000   1.000   1.000  1.00  0.00           N
ATOM      4  CA  HIS A  42       4.000   1.000   1.000  1.00  0.00           C
ATOM      5  N   HIS B  50       5.000   1.000   1.000  1.00  0.00           N
ATOM      6  CA  HIS B  50       6.000   1.000   1.000  1.00  0.00           C
ATOM      7  N   GLU B  51       7.000   1.000   1.000  1.00  0.00           N
ATOM      8  CA  GLU B  51       8.000   1.000   1.000  1.00  0.00           C
END
"""


@pytest.fixture
def mini_pdb(tmp_path):
    p = tmp_path / "test.pdb"
    p.write_text(_MINI_PDB)
    return p


def test_detect_his_residues(mini_pdb):
    his_list = detect_his_residues(mini_pdb)
    assert len(his_list) == 2
    assert his_list[0] == {"chain": "A", "resid": 42, "resname": "HIS"}
    assert his_list[1] == {"chain": "B", "resid": 50, "resname": "HIS"}


def test_filter_titratable_empty_config(mini_pdb):
    his_list = detect_his_residues(mini_pdb)
    filtered = filter_titratable(his_list, [])
    assert len(filtered) == 2  # empty config = keep all


def test_filter_titratable_specific(mini_pdb):
    his_list = detect_his_residues(mini_pdb)
    filtered = filter_titratable(his_list, [{"chain": "A", "resid": 42}])
    assert len(filtered) == 1
    assert filtered[0]["resid"] == 42


def test_build_cphmd_groups():
    his_list = [
        {"chain": "A", "resid": 42, "resname": "HIS"},
        {"chain": "B", "resid": 50, "resname": "HIS"},
    ]
    groups = build_cphmd_groups(his_list, barrier_height=5.0)
    assert len(groups) == 2
    assert groups[0]["name"] == "HIS_A42"
    assert groups[0]["barrier"] == 5.0
    assert groups[1]["name"] == "HIS_B50"
```

- [ ] **Step 3: Run tests**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python -m pytest tests/test_protonation.py -v
```

Expected: All tests PASS.

- [ ] **Step 4: Commit**

```bash
git add third_party/molecular_dynamics/lib/protonation.py \
        third_party/molecular_dynamics/tests/test_protonation.py
git commit -m "feat(md): add His protonation manager for CpHMD setup"
```

---

## Task 5: run_md.py — Execution Entry Point

**Files:**
- Create: `third_party/molecular_dynamics/run_md.py`

- [ ] **Step 1: Write run_md.py**

```python
#!/usr/bin/env python3
# third_party/molecular_dynamics/run_md.py
"""MD 模拟执行入口。

用法：
  # 单 PDB，标准 MD
  python run_md.py --pdb input.pdb

  # 单 PDB，CpHMD 双 pH
  python run_md.py --pdb input.pdb --ph 7.4 6.0

  # 批量执行
  python run_md.py --pdb-dir variants/ --ph 7.4 6.0

  # 指定输出目录
  python run_md.py --pdb input.pdb --ph 6.0 --output-dir experiments/1E62/R3/md/

  # 指定配置文件
  python run_md.py --pdb input.pdb --config configs/my_config.yaml
"""

import argparse
import logging
import sys
from pathlib import Path

# 将模块根目录加入 sys.path
_MODULE_ROOT = Path(__file__).resolve().parent
if str(_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(_MODULE_ROOT))

from lib.config import load_config
from lib.gromacs_wrapper import GromacsWrapper
from lib.protonation import detect_his_residues, filter_titratable, build_cphmd_groups

logger = logging.getLogger(__name__)


def run_single(
    pdb_path: Path,
    ph_values: list[float],
    cfg: dict,
    output_dir: Path,
) -> None:
    """对单个 PDB 运行 MD 模拟。

    如果 ph_values 为空，运行一次标准 MD。
    如果提供 ph_values，为每个 pH 运行一次。
    """
    variant_name = pdb_path.stem
    pdb_path = pdb_path.resolve()

    if not ph_values:
        # 标准 MD，无 pH 指定
        work_dir = output_dir / variant_name
        logger.info("Running standard MD for %s", variant_name)
        wrapper = GromacsWrapper(cfg, work_dir)
        wrapper.run_full_pipeline(pdb_path)
        return

    for ph in ph_values:
        work_dir = output_dir / variant_name / f"pH_{ph}"
        logger.info("Running MD for %s at pH %.1f", variant_name, ph)

        # CpHMD 设置
        cphmd_cfg = cfg["md"].get("cphmd", {})
        if cphmd_cfg.get("enabled", False):
            his_list = detect_his_residues(pdb_path)
            titratable = filter_titratable(
                his_list, cphmd_cfg.get("titratable_residues", [])
            )
            groups = build_cphmd_groups(
                titratable,
                barrier_height=cphmd_cfg["lambda_dynamics"]["barrier_height"],
            )
            logger.info(
                "CpHMD: %d titratable His at pH %.1f", len(titratable), ph
            )
        else:
            groups = []

        wrapper = GromacsWrapper(cfg, work_dir)

        # 如果 CpHMD 启用，需在 wrapper 中设置 groups
        if groups:
            # 将 groups 注入 config 以便 production() 使用
            cfg_copy = _deep_copy_cfg(cfg)
            cfg_copy["md"]["cphmd"]["_runtime_groups"] = groups
            wrapper = GromacsWrapper(cfg_copy, work_dir)

        wrapper.run_full_pipeline(pdb_path, ph=ph)


def _deep_copy_cfg(cfg: dict) -> dict:
    """深拷贝配置（避免修改原始 dict）。"""
    import copy
    return copy.deepcopy(cfg)


def main():
    parser = argparse.ArgumentParser(description="MD 模拟执行")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb", type=Path, help="单个 PDB 文件路径")
    group.add_argument("--pdb-dir", type=Path, help="PDB 文件目录（批量执行）")
    parser.add_argument(
        "--ph", type=float, nargs="+", default=[],
        help="模拟的 pH 值列表（如 7.4 6.0）"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="输出目录（默认从 config 读取）"
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="配置文件路径（默认 configs/md_config.yaml）"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg = load_config(args.config)
    output_dir = args.output_dir or Path(cfg["md"]["output_dir"])
    output_dir = output_dir.resolve()

    # 收集 PDB 文件列表
    if args.pdb:
        pdb_files = [args.pdb]
    else:
        pdb_files = sorted(args.pdb_dir.glob("*.pdb"))
        if not pdb_files:
            logger.error("No PDB files found in %s", args.pdb_dir)
            sys.exit(1)
        logger.info("Found %d PDB files in %s", len(pdb_files), args.pdb_dir)

    for pdb in pdb_files:
        run_single(pdb, args.ph, cfg, output_dir)

    logger.info("All MD simulations complete. Output: %s", output_dir)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Quick smoke test (help message)**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python run_md.py --help
```

Expected: usage message showing `--pdb`, `--ph`, `--output-dir` options.

- [ ] **Step 3: Commit**

```bash
git add third_party/molecular_dynamics/run_md.py
git commit -m "feat(md): add run_md.py execution entry point with batch + CpHMD support"
```

---

## Task 6: Base Analyzer + RMSD Analyzer

**Files:**
- Create: `third_party/molecular_dynamics/lib/analyzers/base.py`
- Create: `third_party/molecular_dynamics/lib/analyzers/rmsd.py`
- Create: `third_party/molecular_dynamics/tests/test_analyzers.py`

First install MDAnalysis, then implement analyzers.

- [ ] **Step 1: Install MDAnalysis**

```bash
conda activate optim-pipe
pip install MDAnalysis
```

- [ ] **Step 2: Write base analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/base.py
"""分析器基类：统一接口 + 收敛判断工具。"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class BaseAnalyzer(ABC):
    """轨迹分析器基类。

    Parameters
    ----------
    traj_dir : Path
        轨迹目录（包含 production.xtc 和 production.tpr）
    cfg : dict
        完整配置
    """

    name: str = "base"

    def __init__(self, traj_dir: Path, cfg: dict):
        self.traj_dir = Path(traj_dir)
        self.cfg = cfg
        self.analysis_cfg = cfg["analysis"]
        self.tpr = self.traj_dir / "production.tpr"
        self.xtc = self.traj_dir / "production.xtc"
        self.output_dir = self.traj_dir / "analysis"
        self.output_dir.mkdir(exist_ok=True)

    def check_inputs(self) -> bool:
        """检查输入文件是否存在。"""
        if not self.tpr.exists():
            logger.error("TPR file not found: %s", self.tpr)
            return False
        if not self.xtc.exists():
            logger.error("XTC file not found: %s", self.xtc)
            return False
        return True

    @abstractmethod
    def run(self) -> dict[str, Any]:
        """执行分析，返回摘要指标。"""

    @abstractmethod
    def save(self, results: dict[str, Any]) -> Path:
        """保存详细结果到 CSV，返回文件路径。"""


def detect_convergence(
    time_ns: np.ndarray,
    values: np.ndarray,
    window_ns: float = 5.0,
    variance_threshold: float = 0.05,
) -> float:
    """检测 RMSD 时间序列的收敛起始时间。

    使用 sliding window 方差判断：当连续窗口的方差 < 阈值时认为收敛。

    Parameters
    ----------
    time_ns : 时间点（ns）
    values : 对应的值（如 RMSD in Å）
    window_ns : 滑动窗口大小（ns）
    variance_threshold : 方差阈值（Å²）

    Returns
    -------
    float  收敛起始时间（ns），未收敛返回 -1.0
    """
    if len(time_ns) < 2:
        return -1.0

    dt = time_ns[1] - time_ns[0]
    window_frames = max(1, int(window_ns / dt))

    for i in range(len(values) - window_frames):
        window = values[i : i + window_frames]
        if np.var(window) < variance_threshold:
            return float(time_ns[i])

    return -1.0
```

- [ ] **Step 3: Write RMSD analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/rmsd.py
"""RMSD 时间序列分析 + 自动收敛检测。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD as MDA_RMSD
import numpy as np
import pandas as pd

from .base import BaseAnalyzer, detect_convergence

logger = logging.getLogger(__name__)


class RMSDAnalyzer(BaseAnalyzer):
    """计算全局 + CDR 区域 RMSD 时间序列。"""

    name = "rmsd"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ca = u.select_atoms("name CA")

        # 全局 RMSD
        rmsd_analysis = MDA_RMSD(ca, ca, ref_frame=0)
        rmsd_analysis.run()
        # results.rmsd: shape (n_frames, 3) — [frame, time(ps), rmsd(Å)]
        time_ns = rmsd_analysis.results.rmsd[:, 1] / 1000.0
        global_rmsd = rmsd_analysis.results.rmsd[:, 2]

        # CDR RMSD
        cdr_rmsds = {}
        cdr_regions = self.analysis_cfg["cdr_regions"]
        for cdr_name, cdr_def in cdr_regions.items():
            sel = (
                f"name CA and segid {cdr_def['chain']} "
                f"and resid {cdr_def['start']}:{cdr_def['end']}"
            )
            # MDAnalysis segid 可能与 chainID 映射不同，尝试 chainID
            cdr_atoms = u.select_atoms(sel)
            if len(cdr_atoms) == 0:
                # 回退到 chainID 选择
                sel = (
                    f"name CA and chainID {cdr_def['chain']} "
                    f"and resid {cdr_def['start']}:{cdr_def['end']}"
                )
                cdr_atoms = u.select_atoms(sel)

            if len(cdr_atoms) == 0:
                logger.warning("CDR %s: no atoms found, skipping", cdr_name)
                cdr_rmsds[cdr_name] = np.full(len(time_ns), np.nan)
                continue

            cdr_rmsd = MDA_RMSD(cdr_atoms, cdr_atoms, ref_frame=0)
            cdr_rmsd.run()
            cdr_rmsds[cdr_name] = cdr_rmsd.results.rmsd[:, 2]

        # 收敛检测
        conv_cfg = self.analysis_cfg["convergence"]
        converge_ns = detect_convergence(
            time_ns, global_rmsd,
            window_ns=conv_cfg["window_ns"],
            variance_threshold=conv_cfg["variance_threshold"],
        )

        # 构建结果
        results = {
            "time_ns": time_ns,
            "global_rmsd": global_rmsd,
            "cdr_rmsds": cdr_rmsds,
            "converge_ns": converge_ns,
            "summary": {
                "rmsd_mean": float(np.mean(global_rmsd)),
                "rmsd_std": float(np.std(global_rmsd)),
                "converge_ns": converge_ns,
            },
        }
        # CDR 摘要
        for name, vals in cdr_rmsds.items():
            key = f"{name.lower()}_rmsd_mean"
            results["summary"][key] = float(np.nanmean(vals))

        self.save(results)
        logger.info(
            "RMSD analysis: mean=%.3f Å, converge=%.1f ns",
            results["summary"]["rmsd_mean"],
            converge_ns,
        )
        return results

    def save(self, results: dict[str, Any]) -> Path:
        """保存 RMSD 时间序列到 CSV。"""
        data = {"time_ns": results["time_ns"], "global": results["global_rmsd"]}
        for name, vals in results["cdr_rmsds"].items():
            data[name] = vals
        df = pd.DataFrame(data)
        out = self.output_dir / "rmsd_timeseries.csv"
        df.to_csv(out, index=False, float_format="%.4f")
        logger.info("Saved: %s", out)
        return out
```

- [ ] **Step 4: Write test for convergence detection**

```python
# third_party/molecular_dynamics/tests/test_analyzers.py
"""Tests for analyzer components."""

import numpy as np
import pytest

from lib.analyzers.base import detect_convergence


def test_detect_convergence_stable():
    """RMSD 一直稳定 → 从 t=0 收敛。"""
    time_ns = np.arange(0, 50, 0.1)  # 50ns, 0.1ns interval
    values = np.ones_like(time_ns) * 2.0 + np.random.normal(0, 0.01, len(time_ns))
    result = detect_convergence(time_ns, values, window_ns=5.0, variance_threshold=0.05)
    assert result >= 0.0
    assert result < 5.0  # should converge early


def test_detect_convergence_drift():
    """RMSD 持续上升 → 不收敛。"""
    time_ns = np.arange(0, 50, 0.1)
    values = time_ns * 0.1  # linearly increasing
    result = detect_convergence(time_ns, values, window_ns=5.0, variance_threshold=0.001)
    assert result == -1.0


def test_detect_convergence_late():
    """RMSD 先漂移后稳定 → 晚期收敛。"""
    time_ns = np.arange(0, 50, 0.1)
    values = np.where(time_ns < 20, time_ns * 0.1, 2.0)
    values += np.random.normal(0, 0.01, len(values))
    result = detect_convergence(time_ns, values, window_ns=5.0, variance_threshold=0.05)
    assert result >= 19.0  # converges around t=20ns
```

- [ ] **Step 5: Run tests**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python -m pytest tests/test_analyzers.py -v
```

Expected: All tests PASS.

- [ ] **Step 6: Update analyzers/__init__.py**

```python
# third_party/molecular_dynamics/lib/analyzers/__init__.py
"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
}
```

- [ ] **Step 7: Commit**

```bash
git add third_party/molecular_dynamics/lib/analyzers/ \
        third_party/molecular_dynamics/tests/test_analyzers.py
git commit -m "feat(md): add base analyzer framework + RMSD analyzer with convergence detection"
```

---

## Task 7: RMSF Analyzer

**Files:**
- Create: `third_party/molecular_dynamics/lib/analyzers/rmsf.py`
- Modify: `third_party/molecular_dynamics/lib/analyzers/__init__.py`

- [ ] **Step 1: Write RMSF analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/rmsf.py
"""Per-residue RMSF（残基柔性）分析。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF as MDA_RMSF
import numpy as np
import pandas as pd

from .base import BaseAnalyzer, detect_convergence

logger = logging.getLogger(__name__)


class RMSFAnalyzer(BaseAnalyzer):
    """计算 per-residue Cα RMSF，重点关注 CDR 和突变位点。"""

    name = "rmsf"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ca = u.select_atoms("name CA")

        # 可选：只取收敛后区间
        # 先快速跑 RMSD 看收敛点
        from MDAnalysis.analysis.rms import RMSD as MDA_RMSD

        rmsd_run = MDA_RMSD(ca, ca, ref_frame=0)
        rmsd_run.run()
        time_ns = rmsd_run.results.rmsd[:, 1] / 1000.0
        global_rmsd = rmsd_run.results.rmsd[:, 2]

        conv_cfg = self.analysis_cfg["convergence"]
        converge_ns = detect_convergence(
            time_ns, global_rmsd,
            window_ns=conv_cfg["window_ns"],
            variance_threshold=conv_cfg["variance_threshold"],
        )

        # 确定分析起始帧
        if converge_ns > 0:
            start_frame = int(np.searchsorted(time_ns, converge_ns))
        else:
            start_frame = 0

        # 对齐到平均结构
        from MDAnalysis.analysis import align
        average = align.AverageStructure(u, ca, ref_frame=0).run(start=start_frame)
        ref = average.results.universe
        align.AlignTraj(u, ref, select="name CA", in_memory=True).run(start=start_frame)

        # 计算 RMSF
        rmsf_analysis = MDA_RMSF(ca).run(start=start_frame)
        rmsf_values = rmsf_analysis.results.rmsf  # Å

        # 收集每个 CA 原子的链、残基号、残基名
        chains = []
        resids = []
        resnames = []
        for atom in ca:
            chains.append(atom.segid if atom.segid else atom.chainID)
            resids.append(atom.resid)
            resnames.append(atom.resname)

        # CDR 区域摘要
        cdr_regions = self.analysis_cfg["cdr_regions"]
        cdr_summary = {}
        for cdr_name, cdr_def in cdr_regions.items():
            mask = [
                (c == cdr_def["chain"] and cdr_def["start"] <= r <= cdr_def["end"])
                for c, r in zip(chains, resids)
            ]
            if any(mask):
                cdr_vals = rmsf_values[mask]
                cdr_summary[cdr_name] = float(np.mean(cdr_vals))
            else:
                cdr_summary[cdr_name] = float("nan")

        results = {
            "chains": chains,
            "resids": resids,
            "resnames": resnames,
            "rmsf": rmsf_values,
            "converge_ns": converge_ns,
            "summary": {
                "rmsf_global_mean": float(np.mean(rmsf_values)),
                "converge_ns": converge_ns,
            },
        }
        for cdr_name, val in cdr_summary.items():
            results["summary"][f"rmsf_{cdr_name.lower()}_mean"] = val

        self.save(results)
        logger.info("RMSF analysis: global mean=%.3f Å", results["summary"]["rmsf_global_mean"])
        return results

    def save(self, results: dict[str, Any]) -> Path:
        df = pd.DataFrame({
            "chain": results["chains"],
            "resid": results["resids"],
            "resname": results["resnames"],
            "rmsf": results["rmsf"],
        })
        out = self.output_dir / "rmsf_per_residue.csv"
        df.to_csv(out, index=False, float_format="%.4f")
        logger.info("Saved: %s", out)
        return out
```

- [ ] **Step 2: Register in __init__.py**

```python
# third_party/molecular_dynamics/lib/analyzers/__init__.py
"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer
from .rmsf import RMSFAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
    "rmsf": RMSFAnalyzer,
}
```

- [ ] **Step 3: Commit**

```bash
git add third_party/molecular_dynamics/lib/analyzers/rmsf.py \
        third_party/molecular_dynamics/lib/analyzers/__init__.py
git commit -m "feat(md): add RMSF per-residue flexibility analyzer"
```

---

## Task 8: Hydrogen Bond Analyzer

**Files:**
- Create: `third_party/molecular_dynamics/lib/analyzers/hbond.py`
- Modify: `third_party/molecular_dynamics/lib/analyzers/__init__.py`

- [ ] **Step 1: Write H-bond analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/hbond.py
"""抗体-抗原界面氢键占有率分析。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


class HBondAnalyzer(BaseAnalyzer):
    """分析抗体-抗原界面的氢键占有率。"""

    name = "hbond"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))

        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]

        # 构建选择语句
        ab_sel = " or ".join(f"segid {c}" for c in ab_chains)
        ag_sel = " or ".join(f"segid {c}" for c in ag_chains)

        # 尝试用 segid，如果选不到原子则用 chainID
        if len(u.select_atoms(f"({ab_sel})")) == 0:
            ab_sel = " or ".join(f"chainID {c}" for c in ab_chains)
            ag_sel = " or ".join(f"chainID {c}" for c in ag_chains)

        # 抗体→抗原 + 抗原→抗体 的氢键
        hb = HydrogenBondAnalysis(
            u,
            donors_sel=f"({ab_sel}) or ({ag_sel})",
            hydrogens_sel=f"({ab_sel}) or ({ag_sel})",
            acceptors_sel=f"({ab_sel}) or ({ag_sel})",
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150,
        )
        hb.run()

        # 只保留跨界面的氢键（供体在 ab、受体在 ag，或反过来）
        results_array = hb.results.hbonds
        if len(results_array) == 0:
            logger.warning("No hydrogen bonds found at interface")
            return self._empty_results()

        # results_array columns: [frame, donor_ix, hydrogen_ix, acceptor_ix, distance, angle]
        ab_indices = set(u.select_atoms(f"({ab_sel})").indices)
        ag_indices = set(u.select_atoms(f"({ag_sel})").indices)

        interface_mask = []
        for row in results_array:
            d_ix = int(row[1])
            a_ix = int(row[3])
            cross = (d_ix in ab_indices and a_ix in ag_indices) or \
                    (d_ix in ag_indices and a_ix in ab_indices)
            interface_mask.append(cross)

        interface_hbonds = results_array[interface_mask]
        n_frames = len(u.trajectory)

        # 占有率：按供体-受体残基对统计
        pair_counts: dict[tuple, int] = {}
        for row in interface_hbonds:
            d_atom = u.atoms[int(row[1])]
            a_atom = u.atoms[int(row[3])]
            pair_key = (
                d_atom.segid or d_atom.chainID, d_atom.resid, d_atom.resname,
                a_atom.segid or a_atom.chainID, a_atom.resid, a_atom.resname,
            )
            frame = int(row[0])
            if pair_key not in pair_counts:
                pair_counts[pair_key] = set()
            pair_counts[pair_key].add(frame)

        # 转为占有率
        occupancy_data = []
        for pair_key, frames in pair_counts.items():
            occupancy_data.append({
                "donor_chain": pair_key[0],
                "donor_resid": pair_key[1],
                "donor_resname": pair_key[2],
                "acceptor_chain": pair_key[3],
                "acceptor_resid": pair_key[4],
                "acceptor_resname": pair_key[5],
                "occupancy": len(frames) / n_frames,
            })
        occupancy_data.sort(key=lambda x: -x["occupancy"])

        # 总界面氢键数时间序列
        frame_counts = np.zeros(n_frames)
        for row in interface_hbonds:
            frame_counts[int(row[0])] += 1

        time_ns = np.array([ts.time / 1000.0 for ts in u.trajectory])

        results = {
            "occupancy": occupancy_data,
            "hbond_timeseries": frame_counts,
            "time_ns": time_ns,
            "summary": {
                "n_hbond_mean": float(np.mean(frame_counts)),
                "n_hbond_std": float(np.std(frame_counts)),
                "n_unique_pairs": len(pair_counts),
            },
        }

        self.save(results)
        logger.info(
            "H-bond analysis: mean=%.1f bonds/frame, %d unique pairs",
            results["summary"]["n_hbond_mean"],
            results["summary"]["n_unique_pairs"],
        )
        return results

    def _empty_results(self) -> dict[str, Any]:
        return {
            "occupancy": [],
            "hbond_timeseries": np.array([]),
            "time_ns": np.array([]),
            "summary": {
                "n_hbond_mean": 0.0,
                "n_hbond_std": 0.0,
                "n_unique_pairs": 0,
            },
        }

    def save(self, results: dict[str, Any]) -> Path:
        # 占有率表
        if results["occupancy"]:
            df = pd.DataFrame(results["occupancy"])
            out = self.output_dir / "hbond_occupancy.csv"
            df.to_csv(out, index=False, float_format="%.4f")
            logger.info("Saved: %s", out)

        # 时间序列
        if len(results["time_ns"]) > 0:
            df_ts = pd.DataFrame({
                "time_ns": results["time_ns"],
                "n_hbonds": results["hbond_timeseries"],
            })
            out_ts = self.output_dir / "hbond_timeseries.csv"
            df_ts.to_csv(out_ts, index=False, float_format="%.4f")

        return self.output_dir / "hbond_occupancy.csv"
```

- [ ] **Step 2: Register in __init__.py**

```python
# third_party/molecular_dynamics/lib/analyzers/__init__.py
"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer
from .rmsf import RMSFAnalyzer
from .hbond import HBondAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
    "rmsf": RMSFAnalyzer,
    "hbond": HBondAnalyzer,
}
```

- [ ] **Step 3: Commit**

```bash
git add third_party/molecular_dynamics/lib/analyzers/hbond.py \
        third_party/molecular_dynamics/lib/analyzers/__init__.py
git commit -m "feat(md): add interface hydrogen bond occupancy analyzer"
```

---

## Task 9: SASA + Contacts Analyzers

**Files:**
- Create: `third_party/molecular_dynamics/lib/analyzers/sasa.py`
- Create: `third_party/molecular_dynamics/lib/analyzers/contacts.py`
- Modify: `third_party/molecular_dynamics/lib/analyzers/__init__.py`

- [ ] **Step 1: Write SASA analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/sasa.py
"""界面 Buried SASA 分析。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


class SASAAnalyzer(BaseAnalyzer):
    """计算 Buried SASA = SASA(ab) + SASA(ag) - SASA(complex)。"""

    name = "sasa"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))

        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]

        ab_sel = " or ".join(f"segid {c}" for c in ab_chains)
        ag_sel = " or ".join(f"segid {c}" for c in ag_chains)

        # 尝试 segid，回退到 chainID
        if len(u.select_atoms(f"({ab_sel})")) == 0:
            ab_sel = " or ".join(f"chainID {c}" for c in ab_chains)
            ag_sel = " or ".join(f"chainID {c}" for c in ag_chains)

        ab_atoms = u.select_atoms(f"({ab_sel})")
        ag_atoms = u.select_atoms(f"({ag_sel})")
        complex_atoms = u.select_atoms(f"({ab_sel}) or ({ag_sel})")

        from MDAnalysis.analysis.sasa import SASA

        # 每 10 帧采样（SASA 计算较慢）
        step = max(1, len(u.trajectory) // 500)

        sasa_ab = SASA(u, ab_atoms).run(step=step)
        sasa_ag = SASA(u, ag_atoms).run(step=step)
        sasa_cx = SASA(u, complex_atoms).run(step=step)

        sasa_ab_vals = sasa_ab.results.sasa
        sasa_ag_vals = sasa_ag.results.sasa
        sasa_cx_vals = sasa_cx.results.sasa
        buried = sasa_ab_vals + sasa_ag_vals - sasa_cx_vals

        time_ns = np.array([
            u.trajectory[i].time / 1000.0
            for i in range(0, len(u.trajectory), step)
        ])

        results = {
            "time_ns": time_ns,
            "sasa_ab": sasa_ab_vals,
            "sasa_ag": sasa_ag_vals,
            "sasa_complex": sasa_cx_vals,
            "buried_sasa": buried,
            "summary": {
                "buried_sasa_mean": float(np.mean(buried)),
                "buried_sasa_std": float(np.std(buried)),
            },
        }

        self.save(results)
        logger.info("SASA analysis: buried mean=%.1f Å²", results["summary"]["buried_sasa_mean"])
        return results

    def save(self, results: dict[str, Any]) -> Path:
        df = pd.DataFrame({
            "time_ns": results["time_ns"],
            "sasa_ab": results["sasa_ab"],
            "sasa_ag": results["sasa_ag"],
            "sasa_complex": results["sasa_complex"],
            "buried_sasa": results["buried_sasa"],
        })
        out = self.output_dir / "sasa_timeseries.csv"
        df.to_csv(out, index=False, float_format="%.2f")
        logger.info("Saved: %s", out)
        return out
```

**注意：** MDAnalysis 的 SASA 模块（`MDAnalysis.analysis.sasa.SASA`）在某些版本中可能名为 `MDAnalysis.analysis.freesasa` 或需要 `freesasa` 包。如果 import 失败，替换为 `from MDAnalysis.analysis import solvation` 或使用 `shrake_rupley` 算法。实现时根据实际安装版本调整 import。

- [ ] **Step 2: Write contacts analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/contacts.py
"""抗体-抗原界面接触数分析。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


class ContactsAnalyzer(BaseAnalyzer):
    """计算抗体-抗原界面的原子接触数（距离 < 阈值）。"""

    name = "contacts"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        cutoff = self.analysis_cfg["contact_distance"]

        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]

        ab_sel = " or ".join(f"segid {c}" for c in ab_chains)
        ag_sel = " or ".join(f"segid {c}" for c in ag_chains)
        if len(u.select_atoms(f"({ab_sel})")) == 0:
            ab_sel = " or ".join(f"chainID {c}" for c in ab_chains)
            ag_sel = " or ".join(f"chainID {c}" for c in ag_chains)

        ab_heavy = u.select_atoms(f"({ab_sel}) and not name H*")
        ag_heavy = u.select_atoms(f"({ag_sel}) and not name H*")

        time_list = []
        contact_counts = []
        # 残基对接触频率
        pair_freq: dict[tuple, int] = {}
        n_frames = 0

        for ts in u.trajectory:
            n_frames += 1
            time_list.append(ts.time / 1000.0)
            dist = distance_array(ab_heavy.positions, ag_heavy.positions)
            contacts = np.argwhere(dist < cutoff)
            contact_counts.append(len(contacts))

            # 残基对统计
            for ab_ix, ag_ix in contacts:
                ab_atom = ab_heavy[ab_ix]
                ag_atom = ag_heavy[ag_ix]
                pair = (
                    ab_atom.segid or ab_atom.chainID, ab_atom.resid,
                    ag_atom.segid or ag_atom.chainID, ag_atom.resid,
                )
                pair_freq[pair] = pair_freq.get(pair, 0) + 1

        # 接触频率矩阵
        contact_matrix = [
            {
                "ab_chain": p[0], "ab_resid": p[1],
                "ag_chain": p[2], "ag_resid": p[3],
                "frequency": count / n_frames,
            }
            for p, count in sorted(pair_freq.items(), key=lambda x: -x[1])
        ]

        results = {
            "time_ns": np.array(time_list),
            "n_contacts": np.array(contact_counts),
            "contact_matrix": contact_matrix,
            "summary": {
                "n_contacts_mean": float(np.mean(contact_counts)),
                "n_contacts_std": float(np.std(contact_counts)),
                "n_unique_pairs": len(pair_freq),
            },
        }

        self.save(results)
        logger.info(
            "Contacts analysis: mean=%.1f contacts/frame",
            results["summary"]["n_contacts_mean"],
        )
        return results

    def save(self, results: dict[str, Any]) -> Path:
        # 时间序列
        df = pd.DataFrame({
            "time_ns": results["time_ns"],
            "n_contacts": results["n_contacts"],
        })
        out = self.output_dir / "contacts_timeseries.csv"
        df.to_csv(out, index=False, float_format="%.4f")

        # 接触频率矩阵
        if results["contact_matrix"]:
            df_matrix = pd.DataFrame(results["contact_matrix"])
            out_matrix = self.output_dir / "contacts_frequency.csv"
            df_matrix.to_csv(out_matrix, index=False, float_format="%.4f")

        logger.info("Saved: %s", out)
        return out
```

- [ ] **Step 3: Update __init__.py**

```python
# third_party/molecular_dynamics/lib/analyzers/__init__.py
"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer
from .rmsf import RMSFAnalyzer
from .hbond import HBondAnalyzer
from .sasa import SASAAnalyzer
from .contacts import ContactsAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
    "rmsf": RMSFAnalyzer,
    "hbond": HBondAnalyzer,
    "sasa": SASAAnalyzer,
    "contacts": ContactsAnalyzer,
}
```

- [ ] **Step 4: Commit**

```bash
git add third_party/molecular_dynamics/lib/analyzers/sasa.py \
        third_party/molecular_dynamics/lib/analyzers/contacts.py \
        third_party/molecular_dynamics/lib/analyzers/__init__.py
git commit -m "feat(md): add SASA buried surface area + interface contacts analyzers"
```

---

## Task 10: MM-GBSA Analyzer (Optional)

**Files:**
- Create: `third_party/molecular_dynamics/lib/analyzers/mmpbsa.py`
- Modify: `third_party/molecular_dynamics/lib/analyzers/__init__.py`

- [ ] **Step 1: Write MM-GBSA analyzer**

```python
# third_party/molecular_dynamics/lib/analyzers/mmpbsa.py
"""MM-PBSA/GBSA 结合自由能计算（可选，需要 gmx_MMPBSA）。"""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .base import BaseAnalyzer, detect_convergence

logger = logging.getLogger(__name__)


class MMPBSAAnalyzer(BaseAnalyzer):
    """使用 gmx_MMPBSA 计算结合自由能。

    需要额外安装：pip install gmx_MMPBSA
    """

    name = "mmpbsa"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        mmpbsa_cfg = self.analysis_cfg.get("mmpbsa", {})
        if not mmpbsa_cfg.get("enabled", False):
            logger.info("MM-PBSA disabled in config, skipping")
            return {"summary": {}}

        # 检查 gmx_MMPBSA 是否可用
        if not shutil.which("gmx_MMPBSA"):
            logger.error("gmx_MMPBSA not found. Install: pip install gmx_MMPBSA")
            return {"summary": {"error": "gmx_MMPBSA not installed"}}

        method = mmpbsa_cfg.get("method", "GBSA")

        # 生成 gmx_MMPBSA 输入文件
        input_file = self._write_mmpbsa_input(method)

        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]

        # 创建 index 文件定义 receptor/ligand 组
        # gmx_MMPBSA 需要 receptor 和 ligand 的 index 组
        idx_file = self._create_index()

        cmd = [
            "gmx_MMPBSA",
            "-i", str(input_file),
            "-cs", str(self.tpr),
            "-ct", str(self.xtc),
            "-ci", str(idx_file),
            "-o", str(self.output_dir / "FINAL_RESULTS_MMPBSA.dat"),
            "-eo", str(self.output_dir / "FINAL_RESULTS_MMPBSA.csv"),
        ]

        logger.info("Running gmx_MMPBSA (%s)...", method)
        result = subprocess.run(
            cmd, cwd=self.output_dir,
            capture_output=True, text=True, timeout=3600 * 4,
        )

        if result.returncode != 0:
            logger.error("gmx_MMPBSA failed:\n%s", result.stderr[-2000:])
            return {"summary": {"error": result.stderr[-500:]}}

        # 解析结果
        return self._parse_results()

    def _write_mmpbsa_input(self, method: str) -> Path:
        """生成 gmx_MMPBSA 输入文件。"""
        frame_interval = self.analysis_cfg.get("mmpbsa", {}).get("frame_interval_ps", 100)

        if method.upper() == "GBSA":
            content = f"""\
&general
  sys_name="mmpbsa",
  startframe=1,
  interval={max(1, frame_interval // 10)},
  forcefields="oldff/leaprc.ff99SB,leaprc.gaff"
/
&gb
  igb=5, saltcon=0.15,
/
"""
        else:
            content = f"""\
&general
  sys_name="mmpbsa",
  startframe=1,
  interval={max(1, frame_interval // 10)},
  forcefields="oldff/leaprc.ff99SB,leaprc.gaff"
/
&pb
  istrng=0.15, fillratio=4.0,
/
"""
        path = self.output_dir / "mmpbsa.in"
        path.write_text(content)
        return path

    def _create_index(self) -> Path:
        """创建 gmx_MMPBSA 所需的 index 文件。"""
        # gmx_MMPBSA 自动从 tpr 中提取 receptor/ligand
        # 如果需要自定义，在这里生成 .ndx 文件
        # 暂时使用自动检测
        path = self.output_dir / "index.ndx"
        # 空 index 让 gmx_MMPBSA 自动生成
        path.write_text("")
        return path

    def _parse_results(self) -> dict[str, Any]:
        """解析 gmx_MMPBSA 输出。"""
        result_file = self.output_dir / "FINAL_RESULTS_MMPBSA.csv"
        if not result_file.exists():
            return {"summary": {"error": "Result file not found"}}

        df = pd.read_csv(result_file)
        dg_values = df["TOTAL"].values if "TOTAL" in df.columns else np.array([])

        results = {
            "dg_values": dg_values,
            "summary": {
                "dG_bind_mean": float(np.mean(dg_values)) if len(dg_values) > 0 else float("nan"),
                "dG_bind_std": float(np.std(dg_values)) if len(dg_values) > 0 else float("nan"),
            },
        }

        self.save(results)
        return results

    def save(self, results: dict[str, Any]) -> Path:
        summary = results.get("summary", {})
        df = pd.DataFrame([summary])
        out = self.output_dir / "mmpbsa_summary.csv"
        df.to_csv(out, index=False, float_format="%.4f")
        logger.info("Saved: %s", out)
        return out
```

- [ ] **Step 2: Update __init__.py**

```python
# third_party/molecular_dynamics/lib/analyzers/__init__.py
"""轨迹分析器注册表。"""

from .rmsd import RMSDAnalyzer
from .rmsf import RMSFAnalyzer
from .hbond import HBondAnalyzer
from .sasa import SASAAnalyzer
from .contacts import ContactsAnalyzer
from .mmpbsa import MMPBSAAnalyzer

ANALYZERS = {
    "rmsd": RMSDAnalyzer,
    "rmsf": RMSFAnalyzer,
    "hbond": HBondAnalyzer,
    "sasa": SASAAnalyzer,
    "contacts": ContactsAnalyzer,
    "mmpbsa": MMPBSAAnalyzer,
}
```

- [ ] **Step 3: Commit**

```bash
git add third_party/molecular_dynamics/lib/analyzers/mmpbsa.py \
        third_party/molecular_dynamics/lib/analyzers/__init__.py
git commit -m "feat(md): add MM-GBSA binding free energy analyzer (optional)"
```

---

## Task 11: analyze_trajectory.py — Analysis Entry Point

**Files:**
- Create: `third_party/molecular_dynamics/analyze_trajectory.py`

- [ ] **Step 1: Write analyze_trajectory.py**

```python
#!/usr/bin/env python3
# third_party/molecular_dynamics/analyze_trajectory.py
"""轨迹分析入口。

用法：
  # 运行所有分析项
  python analyze_trajectory.py --traj experiments/1E62/R3/md/HE1H/pH_7.4/

  # 指定分析项
  python analyze_trajectory.py --traj experiments/1E62/R3/md/HE1H/pH_7.4/ --analyses rmsd rmsf hbond

  # 使用自定义配置
  python analyze_trajectory.py --traj path/to/traj/ --config configs/my_config.yaml
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

_MODULE_ROOT = Path(__file__).resolve().parent
if str(_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(_MODULE_ROOT))

from lib.config import load_config
from lib.analyzers import ANALYZERS

logger = logging.getLogger(__name__)


def run_analyses(
    traj_dir: Path,
    cfg: dict,
    analyses: list[str] | None = None,
) -> dict[str, dict[str, Any]]:
    """运行指定的分析项。

    Parameters
    ----------
    traj_dir : 轨迹目录（含 production.xtc + production.tpr）
    cfg : 完整配置
    analyses : 要运行的分析项名称列表。None = 全部。

    Returns
    -------
    dict  {analyzer_name: results_dict}
    """
    if analyses is None:
        analyses = list(ANALYZERS.keys())

    # mmpbsa 默认不运行除非显式指定
    if "mmpbsa" in analyses and not cfg["analysis"].get("mmpbsa", {}).get("enabled"):
        if analyses == list(ANALYZERS.keys()):
            analyses.remove("mmpbsa")

    all_results = {}
    for name in analyses:
        if name not in ANALYZERS:
            logger.warning("Unknown analyzer: %s (available: %s)", name, list(ANALYZERS.keys()))
            continue
        logger.info("--- Running %s analysis ---", name)
        analyzer = ANALYZERS[name](traj_dir, cfg)
        try:
            results = analyzer.run()
            all_results[name] = results
        except Exception:
            logger.exception("Analyzer %s failed", name)
            all_results[name] = {"summary": {"error": "failed"}}

    # 保存汇总摘要
    summary = {}
    for name, results in all_results.items():
        if "summary" in results:
            for k, v in results["summary"].items():
                summary[f"{name}_{k}"] = v

    summary_path = traj_dir / "analysis" / "summary.json"
    summary_path.parent.mkdir(exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info("Summary saved: %s", summary_path)

    return all_results


def main():
    parser = argparse.ArgumentParser(description="MD 轨迹分析")
    parser.add_argument("--traj", type=Path, required=True, help="轨迹目录")
    parser.add_argument(
        "--analyses", nargs="+", default=None,
        help=f"分析项（可选：{', '.join(ANALYZERS.keys())}）"
    )
    parser.add_argument("--config", type=Path, default=None, help="配置文件路径")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg = load_config(args.config)
    run_analyses(args.traj, cfg, args.analyses)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke test**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python analyze_trajectory.py --help
```

Expected: usage message showing `--traj`, `--analyses`, `--config` options.

- [ ] **Step 3: Commit**

```bash
git add third_party/molecular_dynamics/analyze_trajectory.py
git commit -m "feat(md): add analyze_trajectory.py entry point"
```

---

## Task 12: compare_ph.py — Dual-pH Comparison

**Files:**
- Create: `third_party/molecular_dynamics/compare_ph.py`

- [ ] **Step 1: Write compare_ph.py**

```python
#!/usr/bin/env python3
# third_party/molecular_dynamics/compare_ph.py
"""双 pH 差异比较。

用法：
  python compare_ph.py --variant-dir experiments/1E62/R3/md/HE1H/ --base-ph 7.4 --target-ph 6.0
"""

import argparse
import json
import logging
import sys
from pathlib import Path

_MODULE_ROOT = Path(__file__).resolve().parent
if str(_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(_MODULE_ROOT))

from lib.config import load_config

logger = logging.getLogger(__name__)


def compare_ph(
    variant_dir: Path,
    base_ph: float,
    target_ph: float,
) -> dict:
    """比较同一变体在两个 pH 下的分析结果差异。

    Parameters
    ----------
    variant_dir : 变体目录（含 pH_7.4/ 和 pH_6.0/ 子目录）
    base_ph : 基准 pH（通常 7.4）
    target_ph : 目标 pH（通常 6.0）

    Returns
    -------
    dict  差异指标
    """
    base_dir = variant_dir / f"pH_{base_ph}" / "analysis"
    target_dir = variant_dir / f"pH_{target_ph}" / "analysis"

    if not base_dir.exists():
        raise FileNotFoundError(f"Base pH analysis not found: {base_dir}")
    if not target_dir.exists():
        raise FileNotFoundError(f"Target pH analysis not found: {target_dir}")

    # 加载两个 summary.json
    base_summary = _load_summary(base_dir / "summary.json")
    target_summary = _load_summary(target_dir / "summary.json")

    # 计算差异：Δ = target - base
    deltas = {}
    # 定义要比较的指标
    compare_keys = [
        ("rmsd_rmsd_mean", "delta_rmsd"),
        ("rmsf_rmsf_global_mean", "delta_rmsf_global"),
        ("hbond_n_hbond_mean", "delta_n_hbond"),
        ("sasa_buried_sasa_mean", "delta_buried_sasa"),
        ("contacts_n_contacts_mean", "delta_n_contacts"),
        ("mmpbsa_dG_bind_mean", "delta_dG_bind"),
    ]
    # CDR RMSF
    for cdr in ("h1", "h2", "h3", "l1", "l2", "l3"):
        compare_keys.append(
            (f"rmsf_rmsf_{cdr}_mean", f"delta_rmsf_{cdr}")
        )

    for src_key, delta_key in compare_keys:
        base_val = base_summary.get(src_key)
        target_val = target_summary.get(src_key)
        if base_val is not None and target_val is not None:
            try:
                deltas[delta_key] = float(target_val) - float(base_val)
            except (ValueError, TypeError):
                deltas[delta_key] = None

    result = {
        "variant": variant_dir.name,
        "base_ph": base_ph,
        "target_ph": target_ph,
        "base_summary": base_summary,
        "target_summary": target_summary,
        "deltas": deltas,
    }

    # 保存
    out_path = variant_dir / "ph_comparison.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    logger.info("pH comparison saved: %s", out_path)

    # 打印摘要
    logger.info("=== pH %.1f vs %.1f comparison for %s ===", target_ph, base_ph, variant_dir.name)
    for key, val in deltas.items():
        if val is not None:
            direction = "↑" if val > 0 else "↓" if val < 0 else "="
            logger.info("  %s: %+.4f %s", key, val, direction)

    return result


def _load_summary(path: Path) -> dict:
    if not path.exists():
        logger.warning("Summary not found: %s", path)
        return {}
    with open(path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="双 pH 分析比较")
    parser.add_argument("--variant-dir", type=Path, required=True, help="变体目录")
    parser.add_argument("--base-ph", type=float, default=7.4, help="基准 pH")
    parser.add_argument("--target-ph", type=float, default=6.0, help="目标 pH")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    compare_ph(args.variant_dir, args.base_ph, args.target_ph)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add third_party/molecular_dynamics/compare_ph.py
git commit -m "feat(md): add dual-pH comparison tool"
```

---

## Task 13: report.py — Summary Report Generation

**Files:**
- Create: `third_party/molecular_dynamics/lib/report.py`
- Create: `third_party/molecular_dynamics/tests/test_report.py`

- [ ] **Step 1: Write report.py**

```python
# third_party/molecular_dynamics/lib/report.py
"""汇总所有变体的 MD 分析结果为 md_report.csv。"""

import json
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def generate_report(output_dir: Path, ph_values: list[float] | None = None) -> Path:
    """遍历 output_dir 下所有变体目录，汇总分析结果。

    Parameters
    ----------
    output_dir : MD 输出根目录（如 experiments/1E62/R3/md/）
    ph_values : 要包含的 pH 值列表。None = 自动检测。

    Returns
    -------
    Path  md_report.csv 路径
    """
    output_dir = Path(output_dir)
    rows = []

    for variant_dir in sorted(output_dir.iterdir()):
        if not variant_dir.is_dir():
            continue
        variant_id = variant_dir.name

        # 查找 pH 子目录
        ph_dirs = sorted(variant_dir.glob("pH_*"))
        if not ph_dirs:
            # 可能是标准 MD（无 pH 子目录）
            summary_path = variant_dir / "analysis" / "summary.json"
            if summary_path.exists():
                row = _load_summary_row(summary_path, variant_id, ph=None)
                rows.append(row)
            continue

        for ph_dir in ph_dirs:
            ph_str = ph_dir.name.replace("pH_", "")
            try:
                ph = float(ph_str)
            except ValueError:
                continue
            if ph_values and ph not in ph_values:
                continue

            summary_path = ph_dir / "analysis" / "summary.json"
            if summary_path.exists():
                row = _load_summary_row(summary_path, variant_id, ph=ph)
                rows.append(row)

        # 附加 pH 比较数据
        comparison_path = variant_dir / "ph_comparison.json"
        if comparison_path.exists():
            with open(comparison_path) as f:
                comp = json.load(f)
            deltas = comp.get("deltas", {})
            # 找到对应的行并附加 delta 列
            for row in rows:
                if row["variant_id"] == variant_id and row.get("ph") == comp.get("target_ph"):
                    row.update(deltas)

    if not rows:
        logger.warning("No analysis results found in %s", output_dir)
        return output_dir / "md_report.csv"

    df = pd.DataFrame(rows)
    out = output_dir / "md_report.csv"
    df.to_csv(out, index=False, float_format="%.4f")
    logger.info("Report saved: %s (%d variants)", out, len(df))
    return out


def _load_summary_row(
    summary_path: Path,
    variant_id: str,
    ph: float | None,
) -> dict:
    """从 summary.json 加载一行报告数据。"""
    with open(summary_path) as f:
        summary = json.load(f)
    row = {"variant_id": variant_id, "ph": ph}
    row.update(summary)
    return row
```

- [ ] **Step 2: Write test**

```python
# third_party/molecular_dynamics/tests/test_report.py
"""Tests for report generation."""

import json
from pathlib import Path

import pytest

from lib.report import generate_report


@pytest.fixture
def mock_output_dir(tmp_path):
    """创建模拟的 MD 输出目录结构。"""
    # Variant 1: dual pH
    for ph in ("7.4", "6.0"):
        d = tmp_path / "HE1H" / f"pH_{ph}" / "analysis"
        d.mkdir(parents=True)
        summary = {
            "rmsd_rmsd_mean": 2.0 if ph == "7.4" else 2.5,
            "hbond_n_hbond_mean": 5.0 if ph == "7.4" else 3.0,
        }
        (d / "summary.json").write_text(json.dumps(summary))

    # pH comparison
    comp = {
        "variant": "HE1H",
        "base_ph": 7.4,
        "target_ph": 6.0,
        "deltas": {"delta_rmsd": 0.5, "delta_n_hbond": -2.0},
    }
    (tmp_path / "HE1H" / "ph_comparison.json").write_text(json.dumps(comp))

    # Variant 2: single pH
    d2 = tmp_path / "WT" / "pH_7.4" / "analysis"
    d2.mkdir(parents=True)
    (d2 / "summary.json").write_text(json.dumps({"rmsd_rmsd_mean": 1.5}))

    return tmp_path


def test_generate_report(mock_output_dir):
    report_path = generate_report(mock_output_dir)
    assert report_path.exists()

    import pandas as pd
    df = pd.read_csv(report_path)
    assert len(df) == 3  # HE1H pH7.4, HE1H pH6.0, WT pH7.4
    assert "variant_id" in df.columns
    assert "ph" in df.columns

    # Check delta columns for HE1H pH6.0
    he1h_6 = df[(df["variant_id"] == "HE1H") & (df["ph"] == 6.0)]
    assert len(he1h_6) == 1
    assert he1h_6.iloc[0]["delta_rmsd"] == pytest.approx(0.5)
```

- [ ] **Step 3: Run tests**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python -m pytest tests/test_report.py -v
```

Expected: All tests PASS.

- [ ] **Step 4: Commit**

```bash
git add third_party/molecular_dynamics/lib/report.py \
        third_party/molecular_dynamics/tests/test_report.py
git commit -m "feat(md): add summary report generator"
```

---

## Task 14: Integration Verification

End-to-end verification with real data.

**Prerequisites:** GROMACS must be available, at least one PDB file.

- [ ] **Step 1: Verify GROMACS is accessible**

```bash
# 检查 GROMACS 安装位置和版本
which gmx || module avail gromacs 2>&1 | head -5
gmx --version 2>&1 | head -10
```

If `gmx` is not in PATH, update `md_config.yaml` field `gmx_executable` to the full path. If GROMACS is in a conda environment, note which one and update the config accordingly.

- [ ] **Step 2: Dry-run with a small PDB**

Use a small test PDB (e.g., WT antibody) with very short simulation (1ns production) to verify the full pipeline:

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics

# 临时修改 config 为短模拟
# 将 production.duration_ns 改为 1，便于快速测试

python run_md.py --pdb /path/to/test.pdb --ph 7.4 \
    --output-dir /tmp/md_test/
```

- [ ] **Step 3: Test analysis on the generated trajectory**

```bash
python analyze_trajectory.py \
    --traj /tmp/md_test/<variant>/pH_7.4/ \
    --analyses rmsd rmsf contacts
```

Verify that:
- `analysis/rmsd_timeseries.csv` exists and has data
- `analysis/rmsf_per_residue.csv` exists
- `analysis/contacts_timeseries.csv` exists
- `analysis/summary.json` contains all expected keys

- [ ] **Step 4: Run all unit tests**

```bash
cd /public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics
python -m pytest tests/ -v
```

Expected: All tests PASS.

- [ ] **Step 5: Final commit**

```bash
git add -A third_party/molecular_dynamics/
git commit -m "feat(md): complete MD analysis module — execution automation + 6 analyzers + pH comparison"
```

---

## Verification Plan

1. **Unit tests:** `python -m pytest third_party/molecular_dynamics/tests/ -v` — all pass
2. **Config loading:** verify `md_config.yaml` loads without error
3. **MDP rendering:** verify templates render correctly with config values
4. **GROMACS wrapper:** verify pdb2gmx runs with `-ff`/`-water` name-based selection (no interactive prompts)
5. **Full pipeline:** run a short (1ns) MD simulation on a test PDB
6. **Analysis:** run all analyzers on the generated trajectory and verify CSV outputs
7. **pH comparison:** run dual-pH simulation and verify `compare_ph.py` produces `ph_comparison.json`
8. **Report:** run `generate_report()` and verify `md_report.csv` aggregates all variants
