import matplotlib.pyplot as plt
import numpy as np
import os


def plot_xvg(filepath, basepath, title="XVG Plot", xlabel="X-axis", ylabel="Y-axis", skip_header_char='@#'):
    """
    一个简单的函数来读取和绘制XVG文件。
    通常XVG文件第一列是x，第二列是y。
    """
    x_data = []
    y_data = []

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith(tuple(skip_header_char)): # 跳过注释行和空行
                    try:
                        parts = line.split()
                        x_data.append(float(parts[0]))
                        y_data.append(float(parts[1]))
                    except (ValueError, IndexError):
                        print(f"Skipping malformed line: {line}")
                        continue

        if not x_data or not y_data:
            print(f"No data found in {filepath}")
            return

        plt.figure(figsize=(10, 6))
        plt.plot(x_data, y_data)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.show()
        plt.savefig(os.path.join(basepath, "rmsd_production_ref_em.png"))

    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    basepath = "/public/home/ziyang/code/optim-pipe/data/molecular_dynamics/pdb/fold_2025_04_27_09_07_model_4/"
    filepath = os.path.join(basepath, "rmsd_production_ref_em.xvg")

    plot_xvg(filepath, basepath,
             title="RMSD vs. Time",
             xlabel="Time (ns)",
             ylabel="RMSD (nm)")