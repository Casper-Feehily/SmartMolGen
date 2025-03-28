from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from itertools import combinations
from collections import defaultdict
import re

#/*
#*                        _oo0oo_
#*                       o8888888o
#*                       88" . "88
#*                       (| -_- |)
#*                       0\  =  /0
#*                     ___/`---'\___
#*                   .' \\|     |// '.
#*                  / \\|||  :  |||// \
#*                 / _||||| -:- |||||- \
#*                |   | \\\  - /// |   |
#*                | \_|  ''\---/''  |_/ |
#*                \  .-\__  '-'  ___/-. /
#*              ___'. .'  /--.--\  `. .'___
#*           ."" '<  `.___\_<|>_/___.' >' "".
#*          | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#*          \  \ `_.   \_ __\ /__ _/   .-` /  /
#*      =====`-.____`.___ \_____/___.-`___.-'=====
#*                        `=---='
#*
#*
#*      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*
#*            佛祖保佑       永不宕机     永无BUG

# 计算不饱和度
def calculate_unsaturation(c, h, o=0):
    """计算不饱和度：(2C + 2 - H) / 2"""
    unsat = (2 * c + 2 - h) / 2
    if unsat < 0 or unsat != int(unsat):
        return -1  # 无效分子式
    return int(unsat)


# 生成烷烃异构体
def generate_alkane_isomers(c):
    """生成烷烃骨架的 SMILES"""
    if c <= 0:
        return set()
    if c == 1:
        return {"C"}
    if c == 2:
        return {"CC"}
    if c == 3:
        return {"CCC"}
    if c == 4:
        return {"CCCC", "CC(C)C"}
    if c == 5:
        return {"CCCCC", "CC(C)CC", "C(C)(C)C"}
    # 对于更大的 c，生成链状结构
    return {"C" * c}


# 添加双键或三键
def add_unsaturations(skel_smi, double_bonds, triple_bonds):
    """在骨架上添加指定数量的双键和三键"""
    mol = Chem.MolFromSmiles(skel_smi)
    if not mol:
        return set()
    bonds = [b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE]
    total_unsat = double_bonds + 2 * triple_bonds
    if len(bonds) < double_bonds + triple_bonds:
        return set()
    results = set()
    for bond_combo in combinations(bonds, double_bonds + triple_bonds):
        rw = Chem.RWMol(mol)
        for idx, bond in enumerate(bond_combo[:double_bonds]):
            rw.GetBondWithIdx(bond.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        for idx, bond in enumerate(bond_combo[double_bonds:double_bonds + triple_bonds]):
            rw.GetBondWithIdx(bond.GetIdx()).SetBondType(Chem.BondType.TRIPLE)
        try:
            Chem.SanitizeMol(rw)
            results.add(Chem.MolToSmiles(rw, canonical=True))
        except:
            continue
    return results


# 生成碳骨架
def generate_carbon_skeletons(c, u):
    """生成具有指定碳原子数和不饱和度的骨架"""
    if u < 0 or c <= 0:
        return set()
    skeletons = set()
    # 烷烃骨架
    if u == 0:
        skeletons.update(generate_alkane_isomers(c))
    # 不饱和骨架
    for double in range(u + 1):
        triple = (u - double) // 2
        if double + 2 * triple != u:
            continue
        for skel in generate_alkane_isomers(c):
            unsat_skel = add_unsaturations(skel, double, triple)
            skeletons.update(unsat_skel)
    # 环状骨架
    for ring_size in range(3, min(c + 1, 7)):
        ring_unsat = ring_size - 2  # 单环的不饱和度
        if ring_unsat <= u:
            remaining_c = c - ring_size
            if remaining_c >= 0:
                ring_smi = f"C1{'C' * (ring_size - 1)}1"
                if remaining_c > 0:
                    ring_smi += "C" * remaining_c
                try:
                    mol = Chem.MolFromSmiles(ring_smi)
                    Chem.SanitizeMol(mol)
                    base_smi = Chem.MolToSmiles(mol, canonical=True)
                    skeletons.add(base_smi)
                    # 在环上添加额外不饱和度
                    extra_u = u - ring_unsat
                    if extra_u > 0:
                        for d in range(extra_u + 1):
                            t = (extra_u - d) // 2
                            if d + 2 * t == extra_u:
                                unsat_ring = add_unsaturations(base_smi, d, t)
                                skeletons.update(unsat_ring)
                except:
                    continue
    # 芳香环（苯）
    if c >= 6 and u >= 4:
        benzene = "c1ccccc1" + "C" * (c - 6)
        try:
            mol = Chem.MolFromSmiles(benzene)
            Chem.SanitizeMol(mol)
            skeletons.add(Chem.MolToSmiles(mol, canonical=True))
        except:
            pass
    return skeletons


# 添加官能团通用函数
def add_functional_group(skel_smi, group_type, count):
    """在骨架上添加指定数量的官能团"""
    mol = Chem.MolFromSmiles(skel_smi)
    if not mol:
        return set()
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    if len(carbons) < count:
        return set()
    results = set()
    for positions in combinations(carbons, count):
        rw = Chem.RWMol(mol)
        for pos in positions:
            if group_type == "OH":
                o_idx = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(pos, o_idx, Chem.BondType.SINGLE)
            elif group_type == "CHO":
                c_idx = rw.AddAtom(Chem.Atom(6))
                o_idx = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(pos, c_idx, Chem.BondType.SINGLE)
                rw.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
            elif group_type == "CO":
                o_idx = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(pos, o_idx, Chem.BondType.DOUBLE)
            elif group_type == "COOH":
                c_idx = rw.AddAtom(Chem.Atom(6))
                o1_idx = rw.AddAtom(Chem.Atom(8))
                o2_idx = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(pos, c_idx, Chem.BondType.SINGLE)
                rw.AddBond(c_idx, o1_idx, Chem.BondType.DOUBLE)
                rw.AddBond(c_idx, o2_idx, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
            results.add(Chem.MolToSmiles(rw, canonical=True))
        except:
            continue
    return results


# 生成特定类型化合物
def generate_alcohols(c, h, o):
    """生成醇类"""
    u = calculate_unsaturation(c, h, o)
    if u < 0 or o < 1:
        return set()
    skeletons = generate_carbon_skeletons(c, u)
    alcohols = set()
    for skel in skeletons:
        alcohols.update(add_functional_group(skel, "OH", o))
    return alcohols


def generate_aldehydes(c, h, o):
    """生成醛类"""
    u = calculate_unsaturation(c, h, o)
    if u < 1 or o < 1:
        return set()
    skeletons = generate_carbon_skeletons(c - 1, u - 1)
    aldehydes = set()
    for skel in skeletons:
        aldehydes.update(add_functional_group(skel, "CHO", 1))
    return aldehydes


def generate_ketones(c, h, o):
    """生成酮类"""
    u = calculate_unsaturation(c, h, o)
    if u < 1 or o < 1:
        return set()
    skeletons = generate_carbon_skeletons(c, u - 1)
    ketones = set()
    for skel in skeletons:
        ketones.update(add_functional_group(skel, "CO", 1))
    return ketones


def generate_carboxylic_acids(c, h, o):
    """生成羧酸类"""
    u = calculate_unsaturation(c, h, o)
    if u < 1 or o < 2:
        return set()
    skeletons = generate_carbon_skeletons(c - 1, u - 1)
    acids = set()
    for skel in skeletons:
        acids.update(add_functional_group(skel, "COOH", 1))
    return acids


def generate_esters(c, h, o):
    """生成酯类 R-COO-R'"""
    u = calculate_unsaturation(c, h, o)
    if u < 1 or o < 2:
        return set()
    esters = set()
    for c1 in range(1, c):
        c2 = c - c1 - 1  # -COO- 占用 1 个碳
        if c2 < 0:
            continue
        u1_max = min(u - 1, c1 * 2)  # R 的最大不饱和度
        for u1 in range(u1_max + 1):
            u2 = u - 1 - u1
            if u2 < 0 or u2 > c2 * 2:
                continue
            r1_skel = generate_carbon_skeletons(c1, u1)
            r2_skel = generate_carbon_skeletons(c2, u2)
            for r1 in r1_skel:
                for r2 in r2_skel:
                    smi = f"{r1}C(=O)O{r2}" if r2 else f"{r1}C(=O)O"
                    try:
                        mol = Chem.MolFromSmiles(smi)
                        Chem.SanitizeMol(mol)
                        esters.add(Chem.MolToSmiles(mol, canonical=True))
                    except:
                        continue
    return esters


def generate_ethers(c, h, o):
    """生成醚类 R-O-R'"""
    u = calculate_unsaturation(c, h, o)
    if o < 1:
        return set()
    ethers = set()
    for c1 in range(1, c):
        c2 = c - c1
        u1_max = min(u, c1 * 2)
        for u1 in range(u1_max + 1):
            u2 = u - u1
            if u2 < 0 or u2 > c2 * 2:
                continue
            r1_skel = generate_carbon_skeletons(c1, u1)
            r2_skel = generate_carbon_skeletons(c2, u2)
            for r1 in r1_skel:
                for r2 in r2_skel:
                    smi = f"{r1}O{r2}"
                    try:
                        mol = Chem.MolFromSmiles(smi)
                        Chem.SanitizeMol(mol)
                        ethers.add(Chem.MolToSmiles(mol, canonical=True))
                    except:
                        continue
    return ethers


def generate_phenols(c, h, o):
    """生成酚类（芳香环上的 -OH）"""
    u = calculate_unsaturation(c, h, o)
    if c < 6 or u < 4 or o < 1:
        return set()
    phenols = set()
    base_smi = "c1ccccc1" + "C" * (c - 6)
    try:
        mol = Chem.MolFromSmiles(base_smi)
        Chem.SanitizeMol(mol)
        skel = Chem.MolToSmiles(mol, canonical=True)
        phenols.update(add_functional_group(skel, "OH", o))
    except:
        pass
    return phenols


# 过滤并验证原子数
def filter_by_atom_counts(smiles_set, c, h, o):
    """过滤出符合指定分子式的 SMILES"""
    filtered = set()
    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        formula = rdMolDescriptors.CalcMolFormula(mol)
        match = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        counts = defaultdict(int)
        for atom, num in match:
            counts[atom] = int(num) if num else 1
        if counts["C"] == c and counts["H"] == h and counts["O"] == o:
            filtered.add(smi)
    return filtered


# 分类化合物
def classify_smiles(smiles_set):
    """根据官能团分类化合物"""
    results = defaultdict(set)
    patterns = {
        "烷烃": (
        Chem.MolFromSmarts("[C!H0]"), lambda mol, u: u == 0 and not any(a.GetSymbol() == "O" for a in mol.GetAtoms())),
        "烯烃": (Chem.MolFromSmarts("C=C"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("C=C"))),
        "炔烃": (Chem.MolFromSmarts("C#C"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("C#C"))),
        "芳香化合物": (Chem.MolFromSmarts("c1ccccc1"),
                       lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")) and u >= 4),
        "醇": (Chem.MolFromSmarts("C[OH]"),
               lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("C[OH]")) and not mol.HasSubstructMatch(
                   Chem.MolFromSmarts("c[OH]"))),
        "酚": (Chem.MolFromSmarts("c[OH]"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("c[OH]"))),
        "醚": (Chem.MolFromSmarts("COC"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("COC"))),
        "醛": (Chem.MolFromSmarts("[CH]=O"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("[CH]=O"))),
        "酮": (Chem.MolFromSmarts("CC(=O)C"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)C"))),
        "羧酸": (Chem.MolFromSmarts("C(=O)O"),
                 lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O")) and mol.HasSubstructMatch(
                     Chem.MolFromSmarts("[OH]"))),
        "酯": (Chem.MolFromSmarts("C(=O)OC"), lambda mol, u: mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)OC"))),
    }
    cyclic_patt = Chem.MolFromSmarts("[R]")

    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        u = calculate_unsaturation(*get_atom_counts(smi))
        classified = False
        for name, (patt, cond) in patterns.items():
            if cond(mol, u):
                results[name].add(smi)
                classified = True
                break
        if not classified:
            results["其他"].add(smi)
        if mol.HasSubstructMatch(cyclic_patt) and not any(
                name in ["芳香化合物"] for name in results if smi in results[name]):
            results["环状化合物"].add(smi)
    return results


def get_atom_counts(smi):
    """获取 SMILES 的原子计数"""
    mol = Chem.MolFromSmiles(smi)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    match = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    counts = defaultdict(int)
    for atom, num in match:
        counts[atom] = int(num) if num else 1
    return counts["C"], counts["H"], counts["O"]


# 生成所有化合物
def generate_compounds(c, h, o):
    """生成所有可能的化合物"""
    u = calculate_unsaturation(c, h, o)
    if u < 0:
        return set()
    all_structures = set()
    # 无氧化合物
    if o == 0:
        all_structures.update(generate_carbon_skeletons(c, u))
    # 含氧化合物
    if o >= 1:
        all_structures.update(generate_alcohols(c, h, o))
        all_structures.update(generate_aldehydes(c, h, o))
        all_structures.update(generate_ketones(c, h, o))
        all_structures.update(generate_ethers(c, h, o))
        all_structures.update(generate_phenols(c, h, o))
    if o >= 2:
        all_structures.update(generate_carboxylic_acids(c, h, o))
        all_structures.update(generate_esters(c, h, o))
    return filter_by_atom_counts(all_structures, c, h, o)


# 主函数
def main():
    """程序入口"""
    print("请输入元素数量：")
    try:
        c = int(input("碳原子数 (C): "))
        h = int(input("氢原子数 (H): "))
        o = int(input("氧原子数 (O): "))
        if c > 10 or c < 0 or h < 0 or o < 0:
            print("❌ 输入无效，碳原子数限制为 0-10，所有原子数需非负。")
            return
    except ValueError:
        print("❌ 输入错误，请输入整数。")
        return

    print("\n正在生成结构，请稍候...\n")
    compounds = generate_compounds(c, h, o)
    if not compounds:
        print("❌ 未找到符合条件的结构，可能输入的分子式无效。")
        return
    classified = classify_smiles(compounds)

    print("=== 结构分类结果 ===")
    total = 0
    for category in ["烷烃", "烯烃", "炔烃", "芳香化合物", "醇", "酚", "醚", "醛", "酮", "羧酸", "酯", "环状化合物",
                     "其他"]:
        smiles_list = classified.get(category, set())
        if smiles_list:
            print(f"\n🔹 {category}（{len(smiles_list)} 种）")
            for smi in sorted(smiles_list):
                print(f"  - {smi}")
            total += len(smiles_list)
    print(f"\n✅ 总计结构数：{total}")


if __name__ == "__main__":
    main()