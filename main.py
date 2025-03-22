from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import combinations

# ==============================================================================
# 辅助函数：价键检查与不饱和度计算
# ==============================================================================
def can_add_bond(mol, atom_idx, bond_order):
    """检查指定原子是否能添加指定阶数的键（碳最大4价，氧最大2价）"""
    atom = mol.GetAtomWithIdx(atom_idx)
    current_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
    max_valence = 4 if atom.GetSymbol() == "C" else 2 if atom.GetSymbol() == "O" else 0
    return (current_valence + bond_order) <= max_valence

def calculate_unsaturation(carbon, hydrogen, oxygen):
    """计算不饱和度 (DoU)：(2C + 2 - H - O) / 2"""
    unsat = (2 * carbon + 2 - hydrogen - oxygen) / 2
    if unsat < 0:
        raise ValueError("不饱和度为负，分子式不合理。")
    return unsat

# ==============================================================================
# 碳骨架生成（饱和）
# ==============================================================================
def generate_carbon_skeletons(carbon_count):
    """生成饱和碳骨架（CnH(2n+2)），包括无环和环状"""
    hardcoded = {
        1: {"C"},
        2: {"CC"},
        3: {"CCC"},
        4: {"CCCC", "CC(C)C"},
        5: {"CCCCC", "CC(C)CC", "CC(C)(C)C"},
        6: {"CCCCCC", "CC(C)CCC", "CCC(C)CC", "CC(C)(C)CC", "CC(C)C(C)C"}
    }
    if carbon_count in hardcoded:
        return hardcoded[carbon_count]
    else:
        acyclic = generate_acyclic_skeletons(carbon_count)
        cyclic = generate_cyclic_skeletons(carbon_count)
        return acyclic.union(cyclic)

def generate_acyclic_skeletons(carbon_count, max_depth=10):
    """递归生成无环饱和烷烃骨架"""
    if carbon_count <= 0:
        return set()
    results = set()
    seen = set()
    start = Chem.MolFromSmiles("C")
    seen.add(Chem.MolToSmiles(start, canonical=True))
    _grow_acyclic(start, 1, carbon_count, results, seen, max_depth)
    return results

def _grow_acyclic(mol, current_count, target_count, results, seen, max_depth):
    if current_count == target_count or current_count > max_depth:
        try:
            Chem.SanitizeMol(mol)
            results.add(Chem.MolToSmiles(mol, canonical=True))
        except Exception:
            pass
        return
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() < 4:
            new_mol = Chem.RWMol(mol)
            new_atom_idx = new_mol.AddAtom(Chem.Atom(6))
            new_mol.AddBond(atom.GetIdx(), new_atom_idx, Chem.BondType.SINGLE)
            try:
                Chem.SanitizeMol(new_mol)
                smi = Chem.MolToSmiles(new_mol, canonical=True)
                if smi not in seen:
                    seen.add(smi)
                    _grow_acyclic(new_mol, current_count + 1, target_count, results, seen, max_depth)
            except Exception:
                continue

def generate_cyclic_skeletons(carbon_count):
    """生成简单环状骨架，包括芳香环"""
    results = set()
    if carbon_count >= 3:
        results.add("C1" + "C" * (carbon_count - 1) + "1")  # 单环
    if carbon_count == 4:
        results.add("C1CCC1")  # 环丁烷
    if carbon_count == 5:
        results.add("C1CCCC1")  # 环戊烷
    if carbon_count == 6:
        results.add("c1ccccc1")  # 苯
    if carbon_count == 10:
        results.add("c1ccc2ccccc2c1")  # 萘
    if carbon_count == 9:
        results.add("C1CCC2(C1)CCCC2")  # 螺[4.4]壬烷
    return results

# ==============================================================================
# 不饱和骨架生成（烯烃、炔烃、芳香族）
# ==============================================================================
def generate_unsaturated_variants_extended(skeleton, target_unsat, carbon_count, oxygen_atoms):
    """从骨架生成不饱和变体（双键、三键、环）"""
    results = set()
    try:
        sat_mol = Chem.MolFromSmiles(skeleton)
        if sat_mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(sat_mol)
        Chem.SanitizeMol(sat_mol)
    except Exception:
        return results
    target_H = 2 * carbon_count + 2 - 2 * target_unsat - oxygen_atoms

    # 候选单键（用于双键或三键）
    candidate_bonds = [bond.GetIdx() for bond in sat_mol.GetBonds()
                       if bond.GetBondType() == Chem.BondType.SINGLE and
                       bond.GetBeginAtom().GetSymbol() == "C" and
                       bond.GetEndAtom().GetSymbol() == "C"][:10]

    # 候选环闭合对（拓扑距离 >= 3）
    candidate_ring = []
    dmat = Chem.GetDistanceMatrix(sat_mol)
    for i, j in combinations(range(sat_mol.GetNumAtoms()), 2):
        a_i, a_j = sat_mol.GetAtomWithIdx(i), sat_mol.GetAtomWithIdx(j)
        if a_i.GetSymbol() == "C" and a_j.GetSymbol() == "C" and \
           sat_mol.GetBondBetweenAtoms(i, j) is None and dmat[i][j] >= 3:
            candidate_ring.append((i, j))
    candidate_ring = candidate_ring[:10]

    # 枚举环闭合、双键、三键组合
    for r in range(min(len(candidate_ring) + 1, int(target_unsat) + 1)):
        unsat_needed = target_unsat - r
        ring_combos = list(combinations(candidate_ring, r)) if r > 0 else [()]
        for ring in ring_combos:
            bonds = candidate_bonds
            n = len(bonds)
            for subset_mask in range(1 << n):
                chosen = [bonds[i] for i in range(n) if (subset_mask >> i) & 1]
                for assign_mask in range(1 << len(chosen)):
                    d_count = sum(1 for j in range(len(chosen)) if not (assign_mask >> j) & 1)
                    t_count = sum(1 for j in range(len(chosen)) if (assign_mask >> j) & 1)
                    if d_count + 2 * t_count != unsat_needed:
                        continue
                    D = [chosen[j] for j in range(len(chosen)) if not (assign_mask >> j) & 1]
                    T = [chosen[j] for j in range(len(chosen)) if (assign_mask >> j) & 1]
                    mod_mol = apply_unsat_modifications(sat_mol, D, T, ring)
                    try:
                        Chem.SanitizeMol(mod_mol)
                        molH = Chem.AddHs(mod_mol)
                        h_count = sum(1 for atom in molH.GetAtoms() if atom.GetSymbol() == "H")
                        if h_count == target_H:
                            results.add(Chem.MolToSmiles(mod_mol, canonical=True))
                    except Exception:
                        continue
    return results

def apply_unsat_modifications(mol, double_indices, triple_indices, ring_pairs):
    """应用不饱和修改（双键、三键、环闭合）"""
    mod_mol = Chem.RWMol(mol)
    for bond_idx in double_indices:
        bond = mod_mol.GetBondWithIdx(bond_idx)
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if not (can_add_bond(mod_mol, i, 1) and can_add_bond(mod_mol, j, 1)):
            continue
        mod_mol.RemoveBond(i, j)
        mod_mol.AddBond(i, j, Chem.BondType.DOUBLE)
    for bond_idx in triple_indices:
        bond = mod_mol.GetBondWithIdx(bond_idx)
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if not (can_add_bond(mod_mol, i, 2) and can_add_bond(mod_mol, j, 2)):
            continue
        mod_mol.RemoveBond(i, j)
        mod_mol.AddBond(i, j, Chem.BondType.TRIPLE)
    for i, j in ring_pairs:
        if mod_mol.GetBondBetweenAtoms(i, j) is None and can_add_bond(mod_mol, i, 1) and can_add_bond(mod_mol, j, 1):
            mod_mol.AddBond(i, j, Chem.BondType.SINGLE)
    return mod_mol.GetMol()

# ==============================================================================
# 官能团生成
# ==============================================================================
def add_hydroxyl_groups(skeleton, count=1):
    """在骨架上添加羟基 (-OH)"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ("C", "c")]
    for positions in combinations(indices, count):
        mod_mol = Chem.RWMol(mol)
        valid = True
        for pos in positions:
            if not can_add_bond(mod_mol, pos, 1):
                valid = False
                break
            o_idx = mod_mol.AddAtom(Chem.Atom(8))
            mod_mol.AddBond(pos, o_idx, Chem.BondType.SINGLE)
        if valid:
            try:
                Chem.SanitizeMol(mod_mol)
                results.add(Chem.MolToSmiles(mod_mol))
            except Exception:
                continue
    return results

def generate_aldehydes(skeleton):
    """生成醛类（-CHO）"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetDegree() < 4]
    for idx in indices:
        rw = Chem.RWMol(mol)
        if can_add_bond(rw, idx, 2):
            o_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(idx, o_idx, Chem.BondType.DOUBLE)
            try:
                Chem.SanitizeMol(rw)
                results.add(Chem.MolToSmiles(rw))
            except Exception:
                continue
    return results

def generate_ketones(skeleton):
    """生成酮类（>C=O）"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetDegree() >= 2]
    for idx in indices:
        rw = Chem.RWMol(mol)
        if can_add_bond(rw, idx, 2):
            o_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(idx, o_idx, Chem.BondType.DOUBLE)
            try:
                Chem.SanitizeMol(rw)
                results.add(Chem.MolToSmiles(rw))
            except Exception:
                continue
    return results

def generate_carboxylic_acids(skeleton):
    """生成羧酸（-COOH）"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ("C", "c")]
    for idx in indices:
        rw = Chem.RWMol(mol)
        if can_add_bond(rw, idx, 1):
            new_c = rw.AddAtom(Chem.Atom(6))
            rw.AddBond(idx, new_c, Chem.BondType.SINGLE)
            if can_add_bond(rw, new_c, 2) and can_add_bond(rw, new_c, 1):
                dO = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(new_c, dO, Chem.BondType.DOUBLE)
                sO = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(new_c, sO, Chem.BondType.SINGLE)
                try:
                    Chem.SanitizeMol(rw)
                    results.add(Chem.MolToSmiles(rw))
                except Exception:
                    continue
    return results

def generate_esters(carbon_count):
    """生成酯类（R-COO-R'），R' 至少含一个碳原子"""
    results = set()
    if carbon_count < 3:
        return results
    for i in range(1, carbon_count):
        r1_carbon = i
        r2_carbon = carbon_count - i - 1
        if r2_carbon <= 0:
            continue
        r1_skeletons = generate_carbon_skeletons(r1_carbon)
        r2_skeletons = generate_carbon_skeletons(r2_carbon)
        for r1, r2 in [(r1, r2) for r1 in r1_skeletons for r2 in r2_skeletons]:
            ester_smiles = f"{r1}C(=O)O{r2}"
            try:
                mol = Chem.MolFromSmiles(ester_smiles)
                Chem.SanitizeMol(mol)
                results.add(Chem.MolToSmiles(mol, canonical=True))
            except Exception:
                continue
    return results

def generate_ethers(carbon_count):
    """生成醚类（R-O-R'）"""
    results = set()
    if carbon_count < 2:
        return results
    for i in range(1, carbon_count):
        r1_carbon = i
        r2_carbon = carbon_count - i
        if r2_carbon <= 0:
            continue
        r1_skeletons = generate_carbon_skeletons(r1_carbon)
        r2_skeletons = generate_carbon_skeletons(r2_carbon)
        for r1, r2 in [(r1, r2) for r1 in r1_skeletons for r2 in r2_skeletons]:
            ether_smiles = f"{r1}O{r2}"
            try:
                mol = Chem.MolFromSmiles(ether_smiles)
                Chem.SanitizeMol(mol)
                results.add(Chem.MolToSmiles(mol, canonical=True))
            except Exception:
                continue
    return results

def generate_phenols(aromatic_skeleton, oxygen_count):
    """生成酚类（芳香环上的-OH）"""
    results = set()
    mol = Chem.MolFromSmiles(aromatic_skeleton)
    if not mol or not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return results
    aromatic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic() and atom.GetSymbol() == "C"]
    for positions in combinations(aromatic_carbons, oxygen_count):
        mod_mol = Chem.RWMol(mol)
        valid = True
        for pos in positions:
            if not can_add_bond(mod_mol, pos, 1):
                valid = False
                break
            o_idx = mod_mol.AddAtom(Chem.Atom(8))
            mod_mol.AddBond(pos, o_idx, Chem.BondType.SINGLE)
        if valid:
            try:
                Chem.SanitizeMol(mod_mol)
                results.add(Chem.MolToSmiles(mod_mol))
            except Exception:
                continue
    return results

# ==============================================================================
# 验证、过滤与分类
# ==============================================================================
def validate_smiles(smiles_set):
    """验证 SMILES 有效性，返回 (SMILES, 氢数, 类型)"""
    valid = set()
    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try:
                Chem.SanitizeMol(mol)
                molH = Chem.AddHs(mol)
                h_count = sum(1 for atom in molH.GetAtoms() if atom.GetSymbol() == "H")
                mol_type = classify_compound(smi)
                valid.add((Chem.MolToSmiles(mol, canonical=True), h_count, mol_type))
            except Exception:
                continue
    return valid

def filter_by_hydrogen_count(valid_tuples, target_H):
    """按氢原子数过滤"""
    return [(smi, typ) for smi, h, typ in valid_tuples if h == target_H]

def classify_compound(smiles):
    """使用 SMARTS 模式分类化合物类型"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "未知"
    aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    pat_alcohol = Chem.MolFromSmarts("[CX4]-[OX2H]")
    pat_phenol = Chem.MolFromSmarts("c-[OX2H]")
    pat_aldehyde = Chem.MolFromSmarts("[CX3H1](=O)-[#6]")
    pat_ketone = Chem.MolFromSmarts("[CX3](=O)([#6])[#6]")
    pat_ester = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
    pat_ether = Chem.MolFromSmarts("[#6]-[OX2]-[#6]")
    pat_carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(pat_ester):
        return "酯类"
    elif mol.HasSubstructMatch(pat_carboxyl):
        return "羧酸"
    elif aromatic and mol.HasSubstructMatch(pat_phenol):
        return "酚类"
    elif mol.HasSubstructMatch(pat_alcohol):
        return "醇类"
    elif mol.HasSubstructMatch(pat_aldehyde):
        return "醛类"
    elif mol.HasSubstructMatch(pat_ketone):
        return "酮类"
    elif mol.HasSubstructMatch(pat_ether) and not aromatic:
        return "醚类"
    elif aromatic:
        return "芳香族化合物"
    elif "C#C" in smiles:
        return "炔烃"
    elif "C=C" in smiles or "c=c" in smiles:
        return "烯烃"
    else:
        return "烷烃"

# ==============================================================================
# 根据氧原子数和不饱和度生成化合物
# ==============================================================================
def generate_compounds(base_skeletons, oxygen_count, carbon_count, target_unsat):
    """生成所有可能的化合物"""
    candidates = set()
    aromatic_skeletons = {s for s in base_skeletons if Chem.MolFromSmiles(s) and
                          any(atom.GetIsAromatic() for atom in Chem.MolFromSmiles(s).GetAtoms())}

    if oxygen_count == 0:
        candidates.update(base_skeletons)
    else:
        for skeleton in base_skeletons:
            if oxygen_count >= 1:
                candidates.update(add_hydroxyl_groups(skeleton, min(oxygen_count, 1)))
                candidates.update(generate_aldehydes(skeleton))
                candidates.update(generate_ketones(skeleton))
            if oxygen_count >= 2:
                candidates.update(generate_carboxylic_acids(skeleton))
                candidates.update(add_hydroxyl_groups(skeleton, min(oxygen_count, 2)))
        if oxygen_count == 1:
            candidates.update(generate_ethers(carbon_count))
            for aromatic in aromatic_skeletons:
                candidates.update(generate_phenols(aromatic, 1))
        if oxygen_count == 2:
            candidates.update(generate_esters(carbon_count))
        if oxygen_count > 2:
            for skeleton in base_skeletons:
                candidates.update(add_hydroxyl_groups(skeleton, oxygen_count))
    return candidates

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    try:
        carbon_atoms = int(input("请输入碳原子数: "))
        hydrogen_atoms = int(input("请输入氢原子数: "))
        oxygen_atoms = int(input("请输入氧原子数: "))
    except ValueError:
        print("请输入正确的整数！")
        return

    try:
        target_unsat = calculate_unsaturation(carbon_atoms, hydrogen_atoms, oxygen_atoms)
    except ValueError as e:
        print(e)
        return

    sat_skeletons = generate_carbon_skeletons(carbon_atoms)
    base_skeletons = sat_skeletons if target_unsat == 0 else set()
    if target_unsat > 0:
        for s in sat_skeletons:
            base_skeletons.update(generate_unsaturated_variants_extended(s, target_unsat, carbon_atoms, oxygen_atoms))
        if carbon_atoms == 6:
            base_skeletons.add("c1ccccc1")  # 确保苯环被包含

    # 使用 carbon_atoms 而不是 carbon_count
    compounds = generate_compounds(base_skeletons, oxygen_atoms, carbon_atoms, target_unsat)
    valid = validate_smiles(compounds)
    filtered = filter_by_hydrogen_count(valid, hydrogen_atoms)

    unique = {smi: typ for smi, typ in filtered}

    print(f"\n生成的可能SMILES结构 ({len(unique)} 种):")
    for smi, typ in sorted(unique.items()):
        print(f"{smi} - {typ}")

if __name__ == "__main__":
    main()