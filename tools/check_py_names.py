#!/usr/bin/env python
"""
tools/check_py_names.py — 未定义名字 / 幽灵 import 检查（真正的作用域分析）

**为什么需要它（2026-09-26 实测，台账 E-61）。**
`tools/check_py_syntax.mjs` 原来用一份**写死的名单**（`WATCH_SYMBOLS`）查
「用到 common 的符号却没 import」。名单只有 10 个名字，`common.py` 实际
导出 68 个 —— 于是漏掉 `spot_radius_plot_units` 时本地全绿，CI 跑 25 分钟
才在 `03_spatial_domains.py` 崩掉，`NameError: name 'spot_radius_plot_units'
is not defined`，后续 04–08 步全没跑。

根因不是「名单短了一点」，而是**判据的输入域与它要防的缺陷不匹配**：
漏 import 的名字可以是 common 的任意一个导出。把名单补全成 68 个也只是
把同一种漂移往后推一次（下次 `common.py` 加函数又会漂）。

**为什么现在能做真正的作用域分析。** AGENTS.md 规则 17 当年用写死名单，
理由是「改成 common 导出的所有名字会把函数参数名当用法，要正确处理得做
作用域分析 —— 那是重写一个 linter」。那个顾虑是对的，但作用域分析不必
自己写：标准库 `symtable` 就是 CPython 编译器的符号表，它按作用域给出每个
名字是 local / parameter / imported / free（闭包）/ global，正是这里需要的。

判据：一个名字在某作用域里**被引用**，却既不是该作用域的局部绑定、也不是
闭包自由变量、也没有在模块层绑定、也不是内置 —— 那它**运行时必然 NameError**
（除非有 `import *` / `globals()` / `exec`，这三种本仓库都没有，脚本会先查）。

顺带查**幽灵 import**：`from common import X` 而 common 里没有 X。

**第三条判据（2026-09-26 加，台账 E-70）：循环体里建了图、却没在同一个
循环体里保存。** 这一条不是"名字"问题，而是**缩进错位**问题 —— 而
`python -m py_compile` 只做编译，缩进变化仍然合法，三套图门禁也全都是绿的
（实测：spatial `07_spatial_communication.py` 因为给 `for` 那一行加了缩进、
循环体留在原缩进，3 张图只落盘 1 张，CI + 验收 79 项 + 静态门禁全绿）。
所以判据必须**结构级**（AST），不能扫文本。

用法: python tools/check_py_names.py <repo_root>
退出码 1 表示发现问题。
"""
import ast
import builtins
import symtable
import sys
from pathlib import Path

# 运行时由解释器注入的名字，不算未定义
EXTRA = {
    "__name__", "__file__", "__doc__", "__package__", "__spec__", "__loader__",
    "__builtins__", "__debug__", "WindowsError", "__annotations__", "__dict__",
}

SCAN_DIRS = ("scripts", "tools")

# 逃生舱：文件里确实需要 exec/globals() 时，在**同一行**写
# `# py-names: unsafe-ok —— <一句理由>`。故意做成"要写一句话"的，
# 静默豁免会让检查退化成没有检查（同 tools/check_doc_refs.mjs 的两个逃生舱）。
UNSAFE_OK = "py-names: unsafe-ok"


def iter_py(repo: Path):
    for d in SCAN_DIRS:
        base = repo / d
        if not base.is_dir():
            continue
        for f in sorted(base.rglob("*.py")):
            if "__pycache__" in f.parts:
                continue
            yield f


def scope_first_use(path: Path) -> dict:
    """返回 {(作用域名, 名字): 该作用域内首次引用行号}。

    **必须按作用域定位，不能全文件找首次出现。** 第一版在整棵 AST 上找
    `n.id == nm` 的第一个 `Name`，于是 `plt` 被报成某个**别的函数**里的
    `plt` —— 那个函数自己 `import matplotlib.pyplot as plt`，行号指向一处
    **没问题的代码**。
    指向错的行号比不指行号更糟：下一个人会去读一段正确的代码然后困惑。

    （2026-09-26 更正：这段注释原来带**具体的文件名与行号**。写具体行号时
    它是对的，但行号会随代码改动失效，而失效的引用与指向错的行号是同一类
    毛病 —— 下一个人照着去翻，会翻到别的东西。所以改成不写死行号。
    两仓的这份文件逐字节相同，写死行号在另一仓必然错。）
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out = {}

    def visit(node, scope):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                visit(child, child.name)
            elif isinstance(child, ast.Lambda):
                visit(child, scope)
            else:
                if isinstance(child, ast.Name) and (scope, child.id) not in out:
                    out[(scope, child.id)] = child.lineno
                visit(child, scope)

    visit(tree, "<module>")
    return out


def free_names(path: Path):
    """返回 [(名字, 行号, 作用域说明), ...]，按行号排序。

    用 symtable 做逐作用域分析。**按 (作用域, 名字) 分别记录，不按名字去重** ——
    同一个名字在两个不同作用域里都未定义是两处独立缺陷，去重会漏掉一处。
    """
    src = path.read_text(encoding="utf-8")
    st = symtable.symtable(src, str(path), "exec")
    module_bound = module_level_names_of(path)
    first_use = scope_first_use(path)
    undef = {}

    def walk(table):
        for sym in table.get_symbols():
            nm = sym.get_name()
            if not sym.is_referenced():
                continue
            if sym.is_parameter() or sym.is_imported() or sym.is_assigned():
                continue
            if sym.is_free():
                continue          # 闭包捕获，由外层作用域负责
            if sym.is_declared_global() and nm in module_bound:
                continue
            if nm in module_bound or nm in EXTRA or hasattr(builtins, nm):
                continue
            scope = table.get_name()
            undef[(scope, nm)] = first_use.get((scope, nm), 0)
        for child in table.get_children():
            walk(child)

    walk(st)
    return sorted(((nm, line, scope) for (scope, nm), line in undef.items()),
                  key=lambda t: (t[1], t[2], t[0]))


def module_level_names_of(path: Path) -> set:
    """模块级**所有**绑定名字（含 `import numpy as np` 这类附属绑定）。

    用于**幽灵 import** 检查 —— `from common import np` 在运行期确实能成功
    （common 的命名空间里真有 `np`），所以它不算幽灵。
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                names |= {n.id for n in ast.walk(t) if isinstance(n, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for a in node.names:
                names.add(a.asname or a.name.split(".")[0])
    return names


def module_level_defined_of(path: Path) -> set:
    """模块级**自己定义**的名字（`def` / `class` / 赋值），不含 import 进来的。

    **只用于给出「加进 from common import」的提示。** 第一版用上面那个集合
    给提示，于是删掉 `import numpy as np` 时报出
    「未定义名字 np（common 导出过这个名字 —— 加进 'from common import (...)'）」
    —— `np` 出现在 common 的命名空间里，只是因为 **common 自己也 import 了
    numpy**，它不是 common 提供的东西。照着这条提示改会写出
    `from common import np`，那是把 common 的内部依赖当接口用。

    **提示指向错的地方比不给提示更糟**（同行号归属那处同一个道理）。
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                names |= {n.id for n in ast.walk(t) if isinstance(n, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


# 建图 / 保存调用的尾段名字。用尾段而不是全名，是为了同时认
# `plt.subplots` / `plt.figure` / `fig.savefig` / `plt.close` / `save_fig`。
MAKE_FIG_TAILS = {"subplots", "figure"}
SAVE_FIG_TAILS = {"save_fig", "savefig", "close"}
SC_PL_PREFIX = "sc.pl."


def dotted_name(node) -> str | None:
    """把 `a.b.c(...)` 还原成 `"a.b.c"`；`f(...)` 还原成 `"f"`。"""
    if not isinstance(node, ast.Call):
        return None
    parts = []
    cur = node.func
    while isinstance(cur, ast.Attribute):
        parts.append(cur.attr)
        cur = cur.value
    if isinstance(cur, ast.Name):
        parts.append(cur.id)
        return ".".join(reversed(parts))
    return None


def walk_same_scope(node):
    """遍历子树，但**不进入嵌套的 def / lambda / class** —— 那是另一个作用域，
    在它里面保存 figure 也算"这个循环体保存了"（例如循环里调 `def draw(...)`）。"""
    stack = [node]
    while stack:
        cur = stack.pop()
        yield cur
        for child in ast.iter_child_nodes(cur):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef,
                                  ast.Lambda, ast.ClassDef)):
                continue
            stack.append(child)


def save_helper_names(tree) -> set:
    """**本文件里**「函数体内部含保存调用」的函数名集合。

    没有这一条就会**误报**：`for` 里 `plt.subplots()` 之后调 `_draw(i)`，
    而 `_draw` 内部 `save_fig(...)` —— 循环其实每次都保存了，但
    `walk_same_scope` 不下潜到嵌套 `def`，看不见它。
    **误报会让门禁被关掉**（E-64 首版把 `releases/` 当姊妹仓、E-61 把行号
    指向另一段正确代码，都是同一类"指错地方"）。

    这只是**同文件一层**的近似：跨模块的保存帮助函数仍然看不见
    （那种写法本仓库没有；出现时门禁会判红，把保存挪进循环体即可）。
    """
    out = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for sub in ast.walk(node):
            if isinstance(sub, ast.Call):
                nm = dotted_name(sub)
                if nm and nm.rsplit(".", 1)[-1] in SAVE_FIG_TAILS:
                    out.add(node.name)
                    break
    return out


def figure_loops_without_save(path: Path):
    """返回 `[(建图行号, 循环行号), ...]`：循环体建了图但同一循环体内没有保存。

    **为什么这是一个真缺陷而不是风格问题。** 循环体建图、`save_fig` 掉到循环外时：
      · 循环跑了 N 次，**只有最后一次迭代留下的变量值**被那一次 save 落盘；
      · 于是产物里**少 N-1 张图**，而图名看上去完全合规；
      · 前 N-1 个 figure 没有被关闭（句柄泄漏）。
    `figures:dynamic` 验收项只要求"每组 ≥1 张"（槽位是上限不是精确值），
    `check_figures.mjs` 只查有没有墨 —— **没有一条既有判据看得见它**。
    """
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(path))
    helpers = save_helper_names(tree)
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.For, ast.AsyncFor, ast.While)):
            continue
        makes, saves, mk_line = 0, 0, 0
        for sub in walk_same_scope(node):
            if not isinstance(sub, ast.Call):
                continue
            nm = dotted_name(sub)
            if nm is None:
                continue
            tail = nm.rsplit(".", 1)[-1]
            if tail in MAKE_FIG_TAILS or nm.startswith(SC_PL_PREFIX):
                makes += 1
                mk_line = mk_line or sub.lineno
            if tail in SAVE_FIG_TAILS or tail in helpers:
                saves += 1
        if makes and not saves:
            out.append((mk_line, node.lineno))
    return out


def main() -> int:
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    files = list(iter_py(repo))
    if not files:
        # **没有 Python 文件是"不适用"，不是"失败"。** geo 仓库是纯 R
        # （scripts/ 下 0 个 .py），它跑到这里必须放行 —— 判红会让
        # pre-push 在 geo 上永远红，而那条告警与 geo 的改动毫无关系。
        # 这与"仓库目录找错了"必须长得不一样，所以分开检查。
        if not (repo / "scripts").is_dir() and not (repo / "tools").is_dir():
            print("  既没有 scripts/ 也没有 tools/ —— 仓库根目录传错了？")
            return 1
        print("  本仓库没有 Python 脚本（纯 R 仓库），跳过。")
        return 0

    # 有 import * / globals() / exec 时判据不再可靠，直接拒绝而不是假装通过。
    #
    # **必须用 AST 判定，不能扫子串。** 第一版扫的是裸子串，于是本文件
    # 自己的模式元组 `("import *", "globals()", "exec(")` 命中了它自己 ——
    # 检查器把自己的源码判红。**"检查器要检查的东西"与"检查器描述自己要
    # 检查什么"在文本上无法区分，只有语法结构能区分**（同 AGENTS 规则 25.1：
    # `stripComments` 把字符串换成 `""` 后把要检查的东西本身擦掉了）。
    unsafe = []
    for f in files:
        rel = f.relative_to(repo).as_posix()
        src = f.read_text(encoding="utf-8")
        lines = src.splitlines()
        try:
            tree = ast.parse(src)
        except SyntaxError as e:
            unsafe.append((rel, f"语法错误，无法解析：{e}"))
            continue
        for node in ast.walk(tree):
            why = None
            if isinstance(node, ast.ImportFrom) and any(a.name == "*" for a in node.names):
                why = "import *"
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id in ("globals", "exec", "eval", "vars", "locals"):
                    why = f"{node.func.id}()"
            if why is None:
                continue
            # 逃生舱：同一行写了理由就放过（有意的、且理由可核对）
            ln = lines[node.lineno - 1] if 0 < node.lineno <= len(lines) else ""
            if UNSAFE_OK in ln:
                continue
            unsafe.append((f"{rel}:{node.lineno}", why))
    if unsafe:
        print("  作用域分析在这些文件上不可靠（会让判据失效），请先消除：")
        for loc, why in unsafe:
            print(f"      {loc}: 出现 {why}")
        print(f"  （确实需要时在同一行写 `# {UNSAFE_OK} —— <理由>` 豁免）")
        return 1

    # common 的模块级名字，两个集合分工不同：
    #   common_names   —— 含 import 进来的（`np` 也算）→ 幽灵 import 判据
    #   common_defined —— 只有 common 自己 def/赋值的 → 「加进 import」提示
    common_rel = "scripts/lib/common.py"
    common_path = repo / common_rel
    common_names = module_level_names_of(common_path) if common_path.is_file() else set()
    common_defined = module_level_defined_of(common_path) if common_path.is_file() else set()

    failed = 0
    for f in files:
        rel = f.relative_to(repo).as_posix()
        undef = free_names(f)
        for nm, line, where in undef:
            failed += 1
            loc = f"{rel}:{line}" if line else rel
            hint = ""
            if nm in common_defined:
                hint = "  （common 提供这个名字 —— 加进 'from common import (...)'）"
            print(f"  {loc}  未定义名字 {nm} [{where}]{hint}")

    # 幽灵 import：from common import X 而 common 里没有 X
    for f in files:
        rel = f.relative_to(repo).as_posix()
        tree = ast.parse(f.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "common":
                for a in node.names:
                    nm = a.asname or a.name
                    if common_names and nm not in common_names:
                        failed += 1
                        print(f"  {rel}:{node.lineno}  import 了 common 里不存在的 {nm}")

    # 判据 3（E-70）：循环体建图但不在同一循环体内保存。
    #
    # **为什么不能扫文本**：`for ui, r in ...:` 多缩进 4 个空格、循环体留在
    # 原缩进时，源码里 `save_fig(cfg, f"03-07-02-unit{ui}-{pair_slug}", fig)`
    # 这个**字符串一字不改**，只是归属的块变了。只有 AST 能看见归属。
    for f in files:
        rel = f.relative_to(repo).as_posix()
        for mk_line, loop_line in figure_loops_without_save(f):
            failed += 1
            print(f"  {rel}:{mk_line}  循环体里建了图，但同一个循环体内没有保存"
                  f"（循环在 L{loop_line}）")
            print("      循环跑 N 次只会落盘最后一张，其余 figure 未关闭 —— "
                  "把 save_fig/close 挪进循环体")

    n = len(files)
    if failed:
        print(f"\n{n} 个文件里发现 {failed} 处未定义名字 / 幽灵 import / 循环建图未保存")
        return 1
    print(f"静态结构检查通过（{n} 个文件：逐作用域名字 + 幽灵 import + 循环建图保存）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
