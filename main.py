from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import combinations


# ================================
# 辅助函数：检查是否可以添加键
# ================================
def can_add_bond(mol, atom_idx, bond_order):
    """
    检查指定原子是否有足够空闲的价键以添加指定阶数的键。
    对碳原子采用简单的限制：最大4价（不计算隐式氢）。
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    current_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
    max_valence = 4
    return (current_valence + bond_order) <= max_valence


# ================================
# 生成烷烃异构体的碳骨架
# ================================
# 对于 n<=6，直接使用硬编码的集合；对于 n>6，采用生成树算法（生成所有非同构的烷烃骨架）
def generate_carbon_skeletons(carbon_count):
    """
    生成可能的碳骨架，包括：
      - 烷烃的所有非同构 acyclic（无环）骨架，
      - 以及简单的环状骨架和带支链的环骨架。

    对于 n<=6，采用硬编码；对于 n>6，使用树生成算法来构造所有 acyclic 异构体，
    并用简单规则构造环状骨架（纯环及单一支链环）。

    注意：严格来说，烷烃（CnH2n+2）不包含环状骨架（环烷烃公式为 CnH2n），
    但本函数旨在生成广义上的碳骨架，用于后续官能团引入后过滤。
    """
    # 对于 n<=6，使用预定义集合（硬编码已知的部分异构体）
    alkane_isomers = {
        1: {"C"},
        2: {"CC"},
        3: {"CCC"},
        4: {"CCCC", "CC(C)C"},
        5: {"CCCCC", "CC(C)CC", "CC(C)(C)C"},
        6: {"CCCCCC", "CC(C)CCC", "CCC(C)CC", "CC(C)(C)CC", "CC(C)C(C)C"}
    }
    if carbon_count in alkane_isomers:
        return alkane_isomers[carbon_count]
    else:
        acyclic = generate_acyclic_skeletons(carbon_count)
        cyclic = generate_cyclic_skeletons(carbon_count)
        return acyclic.union(cyclic)


# -------------------------------
# 生成 acyclic（无环）骨架的算法
# -------------------------------
def generate_acyclic_skeletons(carbon_count):
    """
    利用递归生长算法生成所有非同构的无环碳骨架。
    从单个碳原子开始，每次在任一可长成的碳原子上添加一个碳，
    直到达到目标数目。通过 canonical SMILES 去重。
    """
    results = set()
    seen = set()
    start = Chem.MolFromSmiles("C")
    canonical = Chem.MolToSmiles(start, canonical=True)
    seen.add(canonical)
    _grow_acyclic(start, 1, carbon_count, results, seen)
    return results


def _grow_acyclic(mol, current_count, target_count, results, seen):
    if current_count == target_count:
        try:
            Chem.SanitizeMol(mol)
            cano = Chem.MolToSmiles(mol, canonical=True)
            results.add(cano)
        except Exception:
            pass
        return

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() < 4:
            # 对每个可长成的碳原子尝试添加新碳
            new_mol = Chem.RWMol(mol)
            new_atom_idx = new_mol.AddAtom(Chem.Atom(6))
            new_mol.AddBond(atom.GetIdx(), new_atom_idx, Chem.BondType.SINGLE)
            new_mol = new_mol.GetMol()
            try:
                Chem.SanitizeMol(new_mol)
            except Exception:
                continue
            cano = Chem.MolToSmiles(new_mol, canonical=True)
            if cano in seen:
                continue
            seen.add(cano)
            _grow_acyclic(new_mol, current_count + 1, target_count, results, seen)


# -------------------------------
# 生成环状骨架的算法
# -------------------------------
def generate_cyclic_skeletons(carbon_count):
    """
    生成简单的环状骨架，包括：
      1. 纯环：如果 carbon_count >= 3，则生成一个全环（cycloalkane）。
      2. 带单一支链的环：对于 cycle_size 从3到 carbon_count-1，
         计算支链长度 L 使得 cycle_size + L - 1 == carbon_count，
         并生成在环上带一个直链支链的骨架。
    """
    results = set()
    # 1. 纯环：注意纯环的公式为 CnH2n
    if carbon_count >= 3:
        # RDKit 标准 SMILES 表示纯环：例如，7个碳为 "C1CCCCCC1"
        pure_cycle = "C1" + "C" * (carbon_count - 1) + "1"
        results.add(pure_cycle)
    # 2. 带支链的环：cycle_size + L - 1 = carbon_count, L>=2
    for cycle_size in range(3, carbon_count):
        L = carbon_count - cycle_size + 1
        if L < 2:
            continue
        # 生成纯环的 SMILES
        cycle_smiles = "C1" + "C" * (cycle_size - 1) + "1"
        # 将支链固定在环上第一个碳上，用括号表示分支
        branch = "C" * (L - 1)
        # 构造带支链的 SMILES：
        # 例如，对于 cycle_smiles="C1CCCC1" 和 branch="C"，
        # 得到 "C(C)1CCCC1" （注意：此处简单插入分支，可能与环结构连接位置有关）
        modified = "C(" + branch + ")" + cycle_smiles[1:]
        results.add(modified)
    return results


# ================================
# 官能团添加函数（与之前保持一致）
# ================================
def add_hydroxyl_groups(skeleton, count=1):
    """
    在分子骨架中添加羟基(-OH)官能团，用于生成醇或酚类化合物。
    添加前检查目标碳原子是否有足够空闲的价键。
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if mol is None:
        return results
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ("C", "c")]
    for positions in combinations(carbon_indices, count):
        modified_mol = Chem.RWMol(mol)
        valid = True
        for pos in positions:
            if not can_add_bond(modified_mol, pos, 1):
                valid = False
                break
            try:
                o_idx = modified_mol.AddAtom(Chem.Atom(8))
                modified_mol.AddBond(pos, o_idx, Chem.BondType.SINGLE)
            except Exception:
                valid = False
                break
        if not valid:
            continue
        try:
            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            continue
    return results


def add_carbonyl_groups(skeleton, count=1):
    """
    在分子骨架中添加羰基 (C=O) 官能团，用于生成醛或酮化合物。
    添加前检查目标碳原子是否允许形成双键。
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if mol is None:
        return results
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ("C", "c")]
    for positions in combinations(carbon_indices, count):
        modified_mol = Chem.RWMol(mol)
        valid = True
        for pos in positions:
            if not can_add_bond(modified_mol, pos, 2):
                valid = False
                break
            try:
                o_idx = modified_mol.AddAtom(Chem.Atom(8))
                modified_mol.AddBond(pos, o_idx, Chem.BondType.DOUBLE)
            except Exception:
                valid = False
                break
        if not valid:
            continue
        try:
            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            continue
    return results


def generate_carboxylic_acids(skeleton):
    """
    构造羧酸 (-COOH) 化合物，采用两种方法生成，同时检查价键限制。
    """
    results = set()
    if skeleton == "CC":
        results.add("CC(=O)O")
        return results
    mol = Chem.MolFromSmiles(skeleton)
    if mol is None:
        return results
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ("C", "c")]
    for c_idx in carbon_indices:
        try:
            skeleton_mol = Chem.MolFromSmiles(skeleton)
            if skeleton_mol is None:
                continue
            if not can_add_bond(skeleton_mol, c_idx, 1):
                continue
            carboxyl_mol = Chem.MolFromSmiles("C(=O)O")
            if carboxyl_mol is None:
                continue
            combined = Chem.CombineMols(skeleton_mol, carboxyl_mol)
            edit_mol = Chem.RWMol(combined)
            edit_mol.AddBond(c_idx, len(skeleton_mol.GetAtoms()), Chem.BondType.SINGLE)
            final_mol = edit_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            pass
        try:
            modified_mol = Chem.RWMol(mol)
            if not can_add_bond(modified_mol, c_idx, 1):
                continue
            c_idx_carboxyl = modified_mol.AddAtom(Chem.Atom(6))
            modified_mol.AddBond(c_idx, c_idx_carboxyl, Chem.BondType.SINGLE)
            if not can_add_bond(modified_mol, c_idx_carboxyl, 2):
                continue
            o_idx_double = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx_carboxyl, o_idx_double, Chem.BondType.DOUBLE)
            if not can_add_bond(modified_mol, c_idx_carboxyl, 1):
                continue
            o_idx_single = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx_carboxyl, o_idx_single, Chem.BondType.SINGLE)
            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            continue
    return results


def generate_aldehydes(skeleton):
    """
    构造醛类化合物，确保正确添加 CHO 基团，并检查价键条件。
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if mol is None:
        return results
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetDegree() < 4]
    for c_idx in carbon_indices:
        modified_mol = Chem.RWMol(mol)
        if not can_add_bond(modified_mol, c_idx, 2):
            continue
        try:
            o_idx = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            continue
    return results


def generate_ketones(skeleton):
    """
    构造酮类化合物，检查非端基碳原子的价键条件，添加羰基生成酮。
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if mol is None:
        return results
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetDegree() >= 2]
    for c_idx in carbon_indices:
        modified_mol = Chem.RWMol(mol)
        if not can_add_bond(modified_mol, c_idx, 2):
            continue
        try:
            o_idx = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            continue
    return results


def generate_esters(skeleton):
    """
    利用模板法构造酯类化合物，同时检查目标碳原子的价键情况。
    """
    results = set()
    if skeleton == "CC":
        results.add("CC(=O)OC")
        return results
    mol = Chem.MolFromSmiles(skeleton)
    if mol is None:
        return results
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetDegree() < 4]
    for c_idx in carbon_indices:
        modified_mol = Chem.RWMol(mol)
        if not can_add_bond(modified_mol, c_idx, 1):
            continue
        try:
            ester_template = Chem.MolFromSmiles("C(=O)OC")
            if ester_template is None:
                continue
            combined = Chem.CombineMols(modified_mol, ester_template)
            edit_mol = Chem.RWMol(combined)
            edit_mol.AddBond(c_idx, len(modified_mol.GetAtoms()), Chem.BondType.SINGLE)
            final_mol = edit_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except Exception:
            continue
    return results


def generate_ethers(total_carbon):
    """
    根据总碳原子数生成醚类化合物：R-O-R，
    通过将碳链分为两部分构造所有可能组合。
    """
    results = set()
    if total_carbon < 2:
        return results
    for i in range(1, total_carbon):
        group1 = "C" * i
        group2 = "C" * (total_carbon - i)
        ether = group1 + "O" + group2
        results.add(ether)
    return results


# ================================
# 验证与分类函数（基本保持不变）
# ================================
def validate_smiles(smiles_set):
    """
    对每个 SMILES 进行验证，计算显式氢原子数量，生成规范化 SMILES，并分类。
    """
    valid_smiles = set()
    invalid_smiles = set()
    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try:
                Chem.SanitizeMol(mol)
                mol_with_H = Chem.AddHs(mol)
                hydrogen_count = sum(1 for atom in mol_with_H.GetAtoms() if atom.GetSymbol() == "H")
                mol_type = classify_compound(smi)
                mol_no_H = Chem.RemoveHs(mol_with_H)
                canonical_smiles = Chem.MolToSmiles(mol_no_H)
                valid_smiles.add((canonical_smiles, hydrogen_count, mol_type))
            except Exception:
                invalid_smiles.add(smi)
        else:
            invalid_smiles.add(smi)
    return valid_smiles, invalid_smiles


def filter_by_hydrogen_count(valid_smiles_tuples, target_hydrogen_count):
    """
    根据目标氢原子数过滤结构，返回满足条件的 SMILES 与化合物类别。
    """
    filtered_smiles = []
    for smiles, h_count, mol_type in valid_smiles_tuples:
        if h_count == target_hydrogen_count:
            filtered_smiles.append((smiles, mol_type))
    return filtered_smiles


def classify_compound(smiles):
    """
    根据 SMILES 判断化合物类别，利用 SMARTS 模式匹配常见官能团。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "未知"
    is_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if "C(=O)O" in smiles or "OC(=O)" in smiles:
        return "羧酸"
    alcohol_pattern = Chem.MolFromSmarts("[CX4]-[OX2H]")
    phenol_pattern = Chem.MolFromSmarts("c-[OX2H]")
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)-[#6]")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)([#6])[#6]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
    if mol.HasSubstructMatch(phenol_pattern) and is_aromatic:
        return "酚类"
    elif mol.HasSubstructMatch(alcohol_pattern):
        return "醇类"
    elif mol.HasSubstructMatch(aldehyde_pattern):
        return "醛类"
    elif mol.HasSubstructMatch(ketone_pattern):
        return "酮类"
    elif mol.HasSubstructMatch(ester_pattern):
        return "酯类"
    elif is_aromatic:
        return "芳香族化合物"
    else:
        return "烷烃"


def generate_compounds_by_oxygen_count(carbon_skeletons, oxygen_count, total_carbon):
    """
    根据氧原子数量生成化合物：
      - 氧数为0时，直接返回碳骨架（烷烃异构体）；
      - 氧数为1时，生成醇/酚、醛/酮以及醚（基于总碳数）；
      - 氧数为2时，生成羧酸、酯和二元醇；
      - 氧数大于2时，生成多元醇等。
    """
    all_structures = set()
    if oxygen_count == 0:
        return carbon_skeletons
    for skeleton in carbon_skeletons:
        if oxygen_count == 1:
            all_structures.update(add_hydroxyl_groups(skeleton, 1))
            all_structures.update(add_carbonyl_groups(skeleton, 1))
            all_structures.update(generate_ethers(total_carbon))
        elif oxygen_count == 2:
            all_structures.update(generate_carboxylic_acids(skeleton))
            all_structures.update(generate_esters(skeleton))
            all_structures.update(add_hydroxyl_groups(skeleton, 2))
        elif oxygen_count > 2:
            all_structures.update(add_hydroxyl_groups(skeleton, oxygen_count))
    return all_structures


# ================================
# 主函数
# ================================
def main():
    try:
        carbon_atoms = int(input("请输入碳原子数: "))
        hydrogen_atoms = int(input("请输入氢原子数: "))
        oxygen_atoms = int(input("请输入氧原子数: "))
    except ValueError:
        print("请输入正确的整数！")
        return

    # 生成碳骨架：扩展 n>6 时，生成 acyclic 与 cyclic 异构体
    carbon_skeletons = generate_carbon_skeletons(carbon_atoms)
    # 根据氧原子数生成化合物（考虑官能团引入）
    all_structures = generate_compounds_by_oxygen_count(carbon_skeletons, oxygen_atoms, carbon_atoms)
    # 验证 SMILES 并计算氢原子数
    valid_structures, invalid_structures = validate_smiles(all_structures)
    # 根据目标氢原子数过滤
    filtered_structures = filter_by_hydrogen_count(valid_structures, hydrogen_atoms)

    unique_structures = {}
    for smiles, mol_type in filtered_structures:
        unique_structures[smiles] = mol_type

    print(f"\n生成的可能SMILES结构 ({len(unique_structures)} 种):")
    for smiles, mol_type in unique_structures.items():
        print(f"{smiles} - {mol_type}")

    if invalid_structures:
        print(f"\n过滤掉的无效SMILES结构: {len(invalid_structures)}")


if __name__ == "__main__":
    main()