from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import product, combinations


def generate_carbon_skeletons(carbon_count):
    """
    生成不同类型的碳骨架，包括烷烃、烯烃和芳香环
    """
    possible_structures = set()

    # 直链烷烃
    alkane = "C" * carbon_count
    possible_structures.add(alkane)

    # 苯环 (C6H6)
    if carbon_count == 6:
        possible_structures.add("c1ccccc1")

    # 支链烷烃 (对于 C ≥ 4 的情况)
    if carbon_count >= 4:
        # 例如：异丁烷
        branched = "C" * (carbon_count - 3) + "C(C)C"
        possible_structures.add(branched)

    return possible_structures


def add_hydroxyl_groups(skeleton, count=1):
    """
    添加羟基(-OH)官能团
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)

    if mol is None:
        return results

    # 获取碳原子索引
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms()
                      if atom.GetSymbol() in ("C", "c")]

    # 在每个位置添加羟基
    for positions in combinations(carbon_indices, count):
        modified_mol = Chem.RWMol(mol)

        for pos in positions:
            try:
                # 添加氧原子并连接单键
                o_idx = modified_mol.AddAtom(Chem.Atom(8))
                modified_mol.AddBond(pos, o_idx, Chem.BondType.SINGLE)

                try:
                    final_mol = modified_mol.GetMol()
                    smi = Chem.MolToSmiles(final_mol)
                    results.add(smi)
                except:
                    continue
            except:
                continue

    return results


def add_carbonyl_groups(skeleton, count=1):
    """
    添加羰基(C=O)官能团
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)

    if mol is None:
        return results

    # 获取碳原子索引
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms()
                      if atom.GetSymbol() in ("C", "c")]

    # 在每个位置添加羰基
    for positions in combinations(carbon_indices, count):
        modified_mol = Chem.RWMol(mol)

        for pos in positions:
            try:
                # 添加氧原子并连接双键
                o_idx = modified_mol.AddAtom(Chem.Atom(8))
                modified_mol.AddBond(pos, o_idx, Chem.BondType.DOUBLE)

                try:
                    final_mol = modified_mol.GetMol()
                    smi = Chem.MolToSmiles(final_mol)
                    results.add(smi)
                except:
                    continue
            except:
                continue

    return results


def generate_carboxylic_acids(skeleton):
    """
    特殊处理生成羧酸，正确处理-COOH结构
    """
    results = set()

    # 处理特殊情况：C2H4O2 乙酸 (醋酸)
    if skeleton == "CC":
        results.add("CC(=O)O")
        return results

    mol = Chem.MolFromSmiles(skeleton)

    if mol is None:
        return results

    # 获取可能的端基碳原子
    carbon_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ("C", "c"):
            # 对于乙酸这样的小分子，任何碳原子都可以连接羧基
            if len(skeleton) <= 2 or atom.GetDegree() < 4:
                carbon_indices.append(atom.GetIdx())

    for c_idx in carbon_indices:
        try:
            # 方法1：直接构建羧酸SMILES结构
            skeleton_mol = Chem.MolFromSmiles(skeleton)
            if skeleton_mol:
                atom = skeleton_mol.GetAtomWithIdx(c_idx)
                # 复制分子骨架并替换连接的原子
                fragments = Chem.FragmentOnBonds(skeleton_mol, [])
                editable_mol = Chem.EditableMol(fragments)

                # 添加羧基原子
                carboxyl_mol = Chem.MolFromSmiles("C(=O)O")
                if carboxyl_mol:
                    combined = Chem.CombineMols(skeleton_mol, carboxyl_mol)
                    # 创建可编辑分子对象
                    edit_mol = Chem.RWMol(combined)

                    # 连接主链和羧基
                    edit_mol.AddBond(c_idx, len(skeleton_mol.GetAtoms()), Chem.BondType.SINGLE)

                    # 获取最终分子并转换为SMILES
                    final_mol = edit_mol.GetMol()
                    Chem.SanitizeMol(final_mol)
                    smi = Chem.MolToSmiles(final_mol)
                    results.add(smi)
        except Exception as e:
            continue

        # 方法2：分步构建羧酸基团
        try:
            modified_mol = Chem.RWMol(mol)

            # 添加羧基的C=O部分
            c_idx_carboxyl = modified_mol.AddAtom(Chem.Atom(6))  # 添加碳原子
            modified_mol.AddBond(c_idx, c_idx_carboxyl, Chem.BondType.SINGLE)

            # 添加双键氧原子
            o_idx_double = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx_carboxyl, o_idx_double, Chem.BondType.DOUBLE)

            # 添加单键氧原子(羟基)
            o_idx_single = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx_carboxyl, o_idx_single, Chem.BondType.SINGLE)

            # 获取最终分子并转换为SMILES
            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except:
            continue

    return results


def generate_aldehydes(skeleton):
    """
    特殊处理生成醛类化合物，确保正确添加CHO基团
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)

    if mol is None:
        return results

    # 优先选择链端的碳原子
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms()
                      if atom.GetSymbol() == "C" and atom.GetDegree() < 4]

    for c_idx in carbon_indices:
        modified_mol = Chem.RWMol(mol)

        try:
            # 添加醛基 (=O)
            o_idx = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)

            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except:
            continue

    return results


def generate_ketones(skeleton):
    """
    特殊处理生成酮类化合物
    """
    results = set()
    mol = Chem.MolFromSmiles(skeleton)

    if mol is None:
        return results

    # 选择非端基碳原子
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms()
                      if atom.GetSymbol() == "C" and atom.GetDegree() >= 2]

    for c_idx in carbon_indices:
        modified_mol = Chem.RWMol(mol)

        try:
            # 添加酮基 (=O)
            o_idx = modified_mol.AddAtom(Chem.Atom(8))
            modified_mol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)

            final_mol = modified_mol.GetMol()
            Chem.SanitizeMol(final_mol)
            smi = Chem.MolToSmiles(final_mol)
            results.add(smi)
        except:
            continue

    return results


def generate_esters(skeleton):
    """
    生成酯类化合物
    """
    results = set()

    # 处理特殊情况
    if skeleton == "CC":
        results.add("CC(=O)OC")  # 乙酸甲酯
        return results

    mol = Chem.MolFromSmiles(skeleton)

    if mol is None:
        return results

    # 获取可能连接酯基的碳原子
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms()
                      if atom.GetSymbol() == "C" and atom.GetDegree() < 4]

    for c_idx in carbon_indices:
        try:
            # 直接使用模板方法
            ester_template = Chem.MolFromSmiles("C(=O)OC")
            if ester_template and mol:
                combined = Chem.CombineMols(mol, ester_template)
                edit_mol = Chem.RWMol(combined)

                # 连接主链和酯基
                edit_mol.AddBond(c_idx, len(mol.GetAtoms()), Chem.BondType.SINGLE)

                # 获取最终分子
                final_mol = edit_mol.GetMol()
                Chem.SanitizeMol(final_mol)
                smi = Chem.MolToSmiles(final_mol)
                results.add(smi)
        except:
            continue

    return results


def validate_smiles(smiles_set):
    """
    验证SMILES并检查氢原子数量
    """
    valid_smiles = set()
    invalid_smiles = set()

    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            # 计算氢原子数量
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "H")

            # 获取化合物类型
            mol_type = classify_compound(smi)

            # 标准化SMILES (不显示氢原子)
            mol = Chem.RemoveHs(mol)
            canonical_smiles = Chem.MolToSmiles(mol)

            valid_smiles.add((canonical_smiles, hydrogen_count, mol_type))
        else:
            invalid_smiles.add(smi)

    return valid_smiles, invalid_smiles


def filter_by_hydrogen_count(valid_smiles_tuples, target_hydrogen_count):
    """
    根据目标氢原子数量过滤结构
    """
    filtered_smiles = []

    for smiles, h_count, mol_type in valid_smiles_tuples:
        if h_count == target_hydrogen_count:
            filtered_smiles.append((smiles, mol_type))

    return filtered_smiles


def classify_compound(smiles):
    """
    根据SMILES分类化合物
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "未知"

    # 检查是否含有苯环
    is_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())

    # 特别检查羧酸模式
    if "C(=O)O" in smiles or "OC(=O)" in smiles:
        return "羧酸"

    # 使用SMARTS模式匹配官能团
    alcohol_pattern = Chem.MolFromSmarts("[CX4]-[OX2H]")
    phenol_pattern = Chem.MolFromSmarts("c-[OX2H]")
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)-[#6]")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)([#6])[#6]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")

    # 匹配官能团
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return "羧酸"
    elif mol.HasSubstructMatch(phenol_pattern) and is_aromatic:
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


def generate_compounds_by_oxygen_count(carbon_skeletons, oxygen_count):
    """
    根据氧原子数量生成化合物
    """
    all_structures = set()

    if oxygen_count == 0:
        return carbon_skeletons

    for skeleton in carbon_skeletons:
        if oxygen_count == 1:
            # 生成含一个氧原子的化合物
            all_structures.update(add_hydroxyl_groups(skeleton, 1))  # 醇
            all_structures.update(add_carbonyl_groups(skeleton, 1))  # 醛或酮

        elif oxygen_count == 2:
            # 生成含两个氧原子的化合物
            # 羧酸 (-COOH)
            all_structures.update(generate_carboxylic_acids(skeleton))
            # 酯类
            all_structures.update(generate_esters(skeleton))
            # 二元醇
            all_structures.update(add_hydroxyl_groups(skeleton, 2))

        elif oxygen_count > 2:
            # 多元醇等
            all_structures.update(add_hydroxyl_groups(skeleton, oxygen_count))

    return all_structures


def main():
    carbon_atoms = int(input("请输入碳原子数: "))
    hydrogen_atoms = int(input("请输入氢原子数: "))
    oxygen_atoms = int(input("请输入氧原子数: "))

    # 生成碳骨架
    carbon_skeletons = generate_carbon_skeletons(carbon_atoms)

    # 根据氧原子数量生成不同的结构
    all_structures = generate_compounds_by_oxygen_count(carbon_skeletons, oxygen_atoms)

    # 验证SMILES并检查氢原子数量
    valid_structures, invalid_structures = validate_smiles(all_structures)

    # 根据目标氢原子数过滤
    filtered_structures = filter_by_hydrogen_count(valid_structures, hydrogen_atoms)

    # 移除重复结构并排序显示结果
    unique_structures = {}
    for smiles, mol_type in filtered_structures:
        unique_structures[smiles] = mol_type

    # 显示结果
    print(f"\n生成的可能SMILES结构 ({len(unique_structures)} 种):")
    for smiles, mol_type in unique_structures.items():
        print(f"{smiles} - {mol_type}")

    # 显示无效的结构数量
    if invalid_structures:
        print(f"\n过滤掉的无效SMILES结构: {len(invalid_structures)}")


if __name__ == "__main__":
    main()
