from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import combinations, product

# ==============================================================================
# 辅助函数
# ==============================================================================
def can_add_bond(mol, atom_idx, bond_order):
    """检查原子是否可以添加指定键级的键"""
    atom = mol.GetAtomWithIdx(atom_idx)
    current_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
    max_valence = 4 if atom.GetSymbol() == "C" else 2 if atom.GetSymbol() == "O" else 0
    return (current_valence + bond_order) <= max_valence

def calculate_unsaturation(carbon, hydrogen):
    """计算不饱和度 (DoU): (2C + 2 - H) / 2"""
    unsat = (2 * carbon + 2 - hydrogen) / 2
    if unsat < 0 or unsat != int(unsat):
        raise ValueError("无效的分子式（不饱和度为负或非整数）")
    return int(unsat)

def SkeletonDoU(skeleton):
    """计算骨架的不饱和度"""
    mol = Chem.MolFromSmiles(skeleton)
    mol_h = Chem.AddHs(mol)
    h_count = sum(1 for atom in mol_h.GetAtoms() if atom.GetSymbol() == "H")
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
    unsat = (2 * c_count + 2 - h_count) / 2
    return int(unsat)

# ==============================================================================
# 碳骨架生成
# ==============================================================================
def generate_carbon_skeletons(carbon_count):
    """生成饱和碳骨架（无环和环状）"""
    hardcoded = {
        1: {"C"},
        2: {"CC"},
        3: {"CCC", "C(C)C"},
        4: {"CCCC", "CC(C)C", "C(C)(C)C"},
        5: {"CCCCC", "CC(C)CC", "C(C)CCC", "CC(C)(C)C", "C(C)(C)CC"},
        6: {"CCCCCC", "CC(C)CCC", "CCC(C)CC", "CC(C)(C)CC", "CC(C)C(C)C", "C(C)(C)CCC", "c1ccccc1"}  # 包括苯环
    }
    cyclic = {
        3: {"C1CC1"},  # 环丙烷
        4: {"C1CCC1"},  # 环丁烷
        5: {"C1CCCC1"},  # 环戊烷
        6: {"C1CCCCC1"},  # 环己烷
    }
    result = set()
    if carbon_count in hardcoded:
        result.update(hardcoded[carbon_count])
    if carbon_count in cyclic:
        result.update(cyclic[carbon_count])
    return result

# ==============================================================================
# 不饱和变体生成
# ==============================================================================
def generate_unsaturated_variants(skeleton, target_unsat):
    """生成不饱和变体（双键、叁键、环）"""
    results = set([skeleton]) if target_unsat == SkeletonDoU(skeleton) else set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results

    bonds = [b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE and
             b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C"]
    dmat = Chem.GetDistanceMatrix(mol)
    ring_pairs = [(i, j) for i, j in combinations(range(mol.GetNumAtoms()), 2)
                  if mol.GetAtomWithIdx(i).GetSymbol() == "C" and
                  mol.GetAtomWithIdx(j).GetSymbol() == "C" and
                  mol.GetBondBetweenAtoms(i, j) is None and dmat[i][j] >= 3]

    base_unsat = SkeletonDoU(skeleton)
    delta_unsat = target_unsat - base_unsat
    if delta_unsat < 0:
        return results

    def apply_modifications(base_mol, doubles, triples, rings):
        rw = Chem.RWMol(base_mol)
        for b_idx in doubles:
            bond = rw.GetBondWithIdx(b_idx)
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if can_add_bond(rw, i, 1) and can_add_bond(rw, j, 1):
                rw.RemoveBond(i, j)
                rw.AddBond(i, j, Chem.BondType.DOUBLE)
        for b_idx in triples:
            bond = rw.GetBondWithIdx(b_idx)
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if can_add_bond(rw, i, 2) and can_add_bond(rw, j, 2):
                rw.RemoveBond(i, j)
                rw.AddBond(i, j, Chem.BondType.TRIPLE)
        for i, j in rings:
            if rw.GetBondBetweenAtoms(i, j) is None and can_add_bond(rw, i, 1) and can_add_bond(rw, j, 1):
                rw.AddBond(i, j, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw, canonical=True)
        except:
            return None

    for r in range(min(len(ring_pairs) + 1, delta_unsat + 1)):
        remaining_unsat = delta_unsat - r
        for ring_combo in combinations(ring_pairs, r) if r > 0 else [()]:
            for d in range(remaining_unsat + 1):
                t = remaining_unsat - d
                if d + 2 * t == remaining_unsat:
                    for d_bonds in combinations(bonds, d):
                        for t_bonds in combinations([b for b in bonds if b not in d_bonds], t):
                            smi = apply_modifications(mol, [b.GetIdx() for b in d_bonds],
                                                      [b.GetIdx() for b in t_bonds], ring_combo)
                            if smi:
                                results.add(smi)
    return results

# ==============================================================================
# 官能团生成
# ==============================================================================
def add_hydroxyl_groups(skeleton, count):
    """添加羟基 (-OH)，区分醇和酚"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    carbons = [(atom.GetIdx(), atom.GetIsAromatic()) for atom in mol.GetAtoms() if atom.GetSymbol() == "C"]
    for positions in combinations([idx for idx, _ in carbons], count):
        rw = Chem.RWMol(mol)
        valid = True
        for pos in positions:
            if not can_add_bond(rw, pos, 1):
                valid = False
                break
            o_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(pos, o_idx, Chem.BondType.SINGLE)
        if valid:
            try:
                Chem.SanitizeMol(rw)
                results.add(Chem.MolToSmiles(rw, canonical=True))
            except:
                continue
    return results

def generate_aldehydes(skeleton):
    """生成醛 (-CHO)，在骨架末端添加"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and can_add_bond(mol, atom.GetIdx(), 2):
            rw = Chem.RWMol(mol)
            o_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(atom.GetIdx(), o_idx, Chem.BondType.DOUBLE)
            try:
                Chem.SanitizeMol(rw)
                results.add(Chem.MolToSmiles(rw, canonical=True))
            except:
                continue
    return results

def generate_ketones(skeleton):
    """生成酮 (>C=O)，在骨架内部添加"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol:
        return results
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() >= 1 and can_add_bond(mol, atom.GetIdx(), 2):
            rw = Chem.RWMol(mol)
            o_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(atom.GetIdx(), o_idx, Chem.BondType.DOUBLE)
            try:
                Chem.SanitizeMol(rw)
                results.add(Chem.MolToSmiles(rw, canonical=True))
            except:
                continue
    return results

def generate_carboxylic_acids(carbon_count):
    """生成羧酸 (-COOH)，使用 C-1 的骨架"""
    results = set()
    if carbon_count < 1:
        return results
    base_skeletons = generate_carbon_skeletons(carbon_count - 1) if carbon_count > 1 else {""}
    for skeleton in base_skeletons:
        mol = Chem.MolFromSmiles(skeleton) if skeleton else Chem.RWMol()
        if not skeleton:
            c_start = mol.AddAtom(Chem.Atom(6))
        else:
            c_start = 0
        for atom_idx in range(mol.GetNumAtoms()) if skeleton else [c_start]:
            if mol.GetNumAtoms() == 0 or (mol.GetAtomWithIdx(atom_idx).GetSymbol() == "C" and can_add_bond(mol, atom_idx, 1)):
                rw = Chem.RWMol(mol)
                c_idx = rw.AddAtom(Chem.Atom(6))
                if skeleton:
                    rw.AddBond(atom_idx, c_idx, Chem.BondType.SINGLE)
                if can_add_bond(rw, c_idx, 3):
                    o1_idx = rw.AddAtom(Chem.Atom(8))
                    rw.AddBond(c_idx, o1_idx, Chem.BondType.DOUBLE)
                    o2_idx = rw.AddAtom(Chem.Atom(8))
                    rw.AddBond(c_idx, o2_idx, Chem.BondType.SINGLE)
                    try:
                        Chem.SanitizeMol(rw)
                        results.add(Chem.MolToSmiles(rw, canonical=True))
                    except:
                        continue
    return results

def generate_esters(carbon_count):
    """生成酯 (R-COO-R')，调整碳分配"""
    results = set()
    if carbon_count < 2:  # 至少需要2个碳（甲酸甲酯）
        return results
    for r1_c in range(0, carbon_count):  # r1_c=0 表示甲酸部分
        r2_c = carbon_count - r1_c - 1
        if r2_c < 1:
            continue
        r1_skeletons = generate_carbon_skeletons(r1_c) if r1_c > 0 else {""}
        r2_skeletons = generate_carbon_skeletons(r2_c)
        for r1, r2 in product(r1_skeletons, r2_skeletons):
            smi = f"{r1}C(=O)O{r2}" if r1 else f"C(=O)O{r2}"
            try:
                mol = Chem.MolFromSmiles(smi)
                Chem.SanitizeMol(mol)
                results.add(Chem.MolToSmiles(mol, canonical=True))
            except:
                continue
    return results

def generate_ethers(carbon_count):
    """生成醚 (R-O-R')"""
    results = set()
    if carbon_count < 2:
        return results
    for r1_c in range(1, carbon_count):
        r2_c = carbon_count - r1_c
        r1_skeletons = generate_carbon_skeletons(r1_c)
        r2_skeletons = generate_carbon_skeletons(r2_c)
        for r1, r2 in product(r1_skeletons, r2_skeletons):
            smi = f"{r1}O{r2}"
            try:
                mol = Chem.MolFromSmiles(smi)
                Chem.SanitizeMol(mol)
                results.add(Chem.MolToSmiles(mol, canonical=True))
            except:
                continue
    return results

def generate_phenols(skeleton, oxygen_count):
    """生成酚（芳香环上的 -OH）"""
    results = set()
    mol = Chem.MolFromSmiles(skeleton)
    if not mol or not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return results
    aromatic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    for positions in combinations(aromatic_carbons, oxygen_count):
        rw = Chem.RWMol(mol)
        valid = True
        for pos in positions:
            if not can_add_bond(rw, pos, 1):
                valid = False
                break
            o_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(pos, o_idx, Chem.BondType.SINGLE)
        if valid:
            try:
                Chem.SanitizeMol(rw)
                results.add(Chem.MolToSmiles(rw, canonical=True))
            except:
                continue
    return results

# ==============================================================================
# 验证和分类
# ==============================================================================
def validate_and_classify(smiles_set, target_hydrogen, target_oxygen):
    """验证 SMILES 并分类化合物，确保氢和氧原子数正确"""
    results = set()
    patterns = [
        ("羧酸", Chem.MolFromSmarts("[CX3](=O)[OX2H1]")),
        ("酯类", Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")),
        ("酚类", Chem.MolFromSmarts("c-[OX2H]")),
        ("醇类", Chem.MolFromSmarts("[CX4]-[OX2H]")),
        ("醛类", Chem.MolFromSmarts("[CX3H1](=O)")),
        ("酮类", Chem.MolFromSmarts("[CX3](=O)([#6])[#6]")),
        ("醚类", Chem.MolFromSmarts("[#6]-[OX2]-[#6]")),
        ("芳香族化合物", Chem.MolFromSmarts("[c]")),
        ("炔烃", Chem.MolFromSmarts("[C]#[C]")),
        ("烯烃", Chem.MolFromSmarts("[C]=[C]")),
        ("烷烃", Chem.MolFromSmarts("[C!H0]"))
    ]

    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try:
                Chem.SanitizeMol(mol)
                mol_h = Chem.AddHs(mol)
                h_count = sum(1 for atom in mol_h.GetAtoms() if atom.GetSymbol() == "H")
                o_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O")
                if h_count == target_hydrogen and o_count == target_oxygen:
                    typ = "未知"
                    for name, pat in patterns:
                        if mol.HasSubstructMatch(pat):
                            typ = name
                            break
                    results.add((Chem.MolToSmiles(mol, canonical=True), typ))
            except:
                continue
    return results

# ==============================================================================
# 化合物生成
# ==============================================================================
def generate_compounds(carbon_count, oxygen_count, target_unsat, hydrogen_count):
    """生成所有可能的化合物"""
    candidates = set()
    base_skeletons = generate_carbon_skeletons(carbon_count)
    unsaturated_skeletons = set()
    for skeleton in base_skeletons:
        unsaturated_skeletons.update(generate_unsaturated_variants(skeleton, target_unsat))

    # 无氧：仅生成烃类
    if oxygen_count == 0:
        candidates.update(unsaturated_skeletons)
    else:
        aromatic_skeletons = {s for s in unsaturated_skeletons if
                              any(atom.GetIsAromatic() for atom in Chem.MolFromSmiles(s).GetAtoms())}

        # 氧 = 1
        if oxygen_count >= 1:
            for skeleton in unsaturated_skeletons:
                candidates.update(add_hydroxyl_groups(skeleton, 1))
                candidates.update(generate_aldehydes(skeleton))
                candidates.update(generate_ketones(skeleton))
            candidates.update(generate_ethers(carbon_count))
            for aromatic in aromatic_skeletons:
                candidates.update(generate_phenols(aromatic, 1))

        # 氧 = 2
        if oxygen_count >= 2:
            candidates.update(generate_carboxylic_acids(carbon_count))
            candidates.update(generate_esters(carbon_count))
            for skeleton in unsaturated_skeletons:
                candidates.update(add_hydroxyl_groups(skeleton, 2))

        # 氧 > 2
        if oxygen_count > 2:
            for skeleton in unsaturated_skeletons:
                candidates.update(add_hydroxyl_groups(skeleton, oxygen_count))

    return validate_and_classify(candidates, hydrogen_count, oxygen_count)

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
        target_unsat = calculate_unsaturation(carbon_atoms, hydrogen_atoms)
    except ValueError as e:
        print(e)
        return

    compounds = generate_compounds(carbon_atoms, oxygen_atoms, target_unsat, hydrogen_atoms)
    print(f"\n生成的可能SMILES结构 ({len(compounds)} 种):")
    for smi, typ in sorted(compounds):
        print(f"{smi} - {typ}")

if __name__ == "__main__":
    main()