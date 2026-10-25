"""DeepTernary 推理子任务（在 deepternary 独立环境内运行）。

用法：deepternary_python dt_task.py <task.json>
task.json: {"repo": ..., "tasks": [{"name", "input_dir", "out_dir", "seeds"}], "result_path": ...}
结果写 result_path（JSON）：每项 {name, summary_csv | error}。
模型只加载一次，批量推理所有分子。
"""
import json
import os
import sys


def main():
    task_path = sys.argv[1]
    task = json.load(open(task_path, encoding="utf-8"))
    repo = task["repo"]
    os.chdir(repo)
    sys.path.insert(0, repo)

    import predict as P  # noqa: E402  cwd=repo 以加载 vendored 包

    cfg, model = P._load_model("deepternary/configs/protac.py",
                               "output/checkpoints/PROTAC/checkpoint.pth", "cpu")
    results = []
    for t in task["tasks"]:
        try:
            summary = P._run_protac(
                name=t["name"],
                lig_path=os.path.join(t["input_dir"], "lig.pdb"),
                p1_path=os.path.join(t["input_dir"], "p1.pdb"),
                p2_path=os.path.join(t["input_dir"], "p2.pdb"),
                lig1_mask_path=os.path.join(t["input_dir"], "lig1_mask.pdb"),
                lig2_mask_path=os.path.join(t["input_dir"], "lig2_mask.pdb"),
                unbound_lig1_path=os.path.join(t["input_dir"], "unbound_lig1.pdb"),
                unbound_lig2_path=os.path.join(t["input_dir"], "unbound_lig2.pdb"),
                outdir=t["out_dir"],
                cfg=cfg,
                num_seeds=t.get("seeds", 2),
                apply_correction=True,
            )
            results.append({"name": t["name"], "summary_csv": summary})
            print(f"[ternary-task] {t['name']} ok -> {summary}", flush=True)
        except Exception as e:  # noqa: BLE001 —— 单分子失败不阻塞批次
            results.append({"name": t["name"], "error": f"{type(e).__name__}: {e}"[:300]})
            print(f"[ternary-task] {t['name']} FAILED: {e}", flush=True)

    with open(task["result_path"], "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=1)
    print("[ternary-task] done", flush=True)


if __name__ == "__main__":
    main()
