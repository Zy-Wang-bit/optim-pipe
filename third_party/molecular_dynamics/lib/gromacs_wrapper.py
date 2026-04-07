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
